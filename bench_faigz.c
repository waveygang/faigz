#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>

#define REENTRANT_FAIDX_IMPLEMENTATION
#include "faigz.h"

// Benchmark configuration
typedef struct {
    char *fasta_file;         // Path to input FASTA file
    int num_threads;          // Number of threads to use
    int seq_count;            // Number of sequences to fetch per thread
    int seq_length;           // Length of sequences to fetch
    char *output_file;        // Optional output file (NULL for no output)
    unsigned int seed;        // PRNG seed
    int verbose;              // Verbose output
} bench_config_t;

// Per-thread data
typedef struct {
    int thread_id;                // Thread ID
    faidx_meta_t *meta;           // Shared metadata
    const bench_config_t *config; // Benchmark configuration
    uint64_t num_bases;           // Total bases retrieved
    double elapsed_time;          // Time spent in seconds
    unsigned int seed;            // Thread-specific PRNG seed
    pthread_mutex_t *output_mutex; // Mutex for writing to output
    FILE *output_fp;              // Output file (shared)
} thread_data_t;

// Function to display usage info
void usage(const char *prog) {
    fprintf(stderr, 
        "Usage: %s [options] <fasta_file>\n"
        "Options:\n"
        "  -t INT    Number of threads [4]\n"
        "  -n INT    Number of sequences to fetch per thread [1000]\n"
        "  -l INT    Length of each sequence to fetch [100]\n"
        "  -o FILE   Output fetched sequences to file [none]\n"
        "  -s INT    Random seed [42]\n"
        "  -v        Verbose output\n"
        "  -h        Show this help message\n",
        prog
    );
}

// Parse command line arguments
bench_config_t parse_args(int argc, char **argv) {
    bench_config_t config = {
        .fasta_file = NULL,
        .num_threads = 4,
        .seq_count = 1000,
        .seq_length = 100,
        .output_file = NULL,
        .seed = 42,
        .verbose = 0
    };

    int c;
    while ((c = getopt(argc, argv, "t:n:l:o:s:vh")) != -1) {
        switch (c) {
            case 't': config.num_threads = atoi(optarg); break;
            case 'n': config.seq_count = atoi(optarg); break;
            case 'l': config.seq_length = atoi(optarg); break;
            case 'o': config.output_file = optarg; break;
            case 's': config.seed = atoi(optarg); break;
            case 'v': config.verbose = 1; break;
            case 'h': usage(argv[0]); exit(0);
            default: usage(argv[0]); exit(1);
        }
    }

    if (optind >= argc) {
        fprintf(stderr, "Error: No FASTA file specified\n");
        usage(argv[0]);
        exit(1);
    }

    config.fasta_file = argv[optind];

    // Validate parameters
    if (config.num_threads < 1) {
        fprintf(stderr, "Error: Number of threads must be >= 1\n");
        exit(1);
    }
    if (config.seq_count < 1) {
        fprintf(stderr, "Error: Number of sequences must be >= 1\n");
        exit(1);
    }
    if (config.seq_length < 1) {
        fprintf(stderr, "Error: Sequence length must be >= 1\n");
        exit(1);
    }

    return config;
}

// Thread worker function
void *worker_thread(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    faidx_reader_t *reader;
    char *seq;
    hts_pos_t seq_len;
    struct timespec start_time, end_time;
    int i, seq_idx;
    uint64_t bases_fetched = 0;
    
    // Initialize thread-specific PRNG
    unsigned int thread_seed = data->seed + data->thread_id;
    srand(thread_seed);
    
    if (data->config->verbose) {
        printf("Thread %d: Starting with seed %u\n", data->thread_id, thread_seed);
    }
    
    // Create a reader
    reader = faidx_reader_create(data->meta);
    if (!reader) {
        fprintf(stderr, "Thread %d: Failed to create reader\n", data->thread_id);
        pthread_exit(NULL);
    }
    
    // Get the number of sequences in the file
    int num_seqs = faidx_meta_nseq(data->meta);
    if (num_seqs <= 0) {
        fprintf(stderr, "Thread %d: No sequences found in the file\n", data->thread_id);
        faidx_reader_destroy(reader);
        pthread_exit(NULL);
    }
    
    // Start timing
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    
    // Fetch random sequences
    for (i = 0; i < data->config->seq_count; i++) {
        // Choose a random sequence
        seq_idx = rand() % num_seqs;
        const char *seq_name = faidx_meta_iseq(data->meta, seq_idx);
        
        // Get the sequence length
        hts_pos_t total_seq_len = faidx_meta_seq_len(data->meta, seq_name);
        if (total_seq_len <= 0) {
            if (data->config->verbose) {
                fprintf(stderr, "Thread %d: Invalid sequence length for %s, skipping\n", 
                       data->thread_id, seq_name);
            }
            continue;
        }
        
        // Choose a random start position, ensuring we don't go past the end
        // Adjust sequence length if the sequence is shorter than requested
        int adjusted_seq_length = data->config->seq_length;
        if (total_seq_len < adjusted_seq_length) {
            adjusted_seq_length = total_seq_len;
            if (data->config->verbose) {
                printf("Thread %d: Sequence %s is shorter than requested length (%"PRIhts_pos" < %d), adjusting\n",
                      data->thread_id, seq_name, total_seq_len, data->config->seq_length);
            }
        }
        
        hts_pos_t max_start = total_seq_len - adjusted_seq_length;
        if (max_start < 0) max_start = 0;
        
        hts_pos_t start = max_start > 0 ? rand() % (max_start + 1) : 0;
        hts_pos_t end = start + adjusted_seq_length - 1;
        if (end >= total_seq_len) end = total_seq_len - 1;
        
        // Fetch the sequence
        seq = faidx_reader_fetch_seq(reader, seq_name, start, end, &seq_len);
        if (!seq) {
            if (data->config->verbose) {
                fprintf(stderr, "Thread %d: Failed to fetch %s:%"PRIhts_pos"-%"PRIhts_pos"\n", 
                       data->thread_id, seq_name, start, end);
            }
            continue;
        }
        
        // Update the number of bases fetched
        bases_fetched += seq_len;
        
        // Write to output file if requested
        if (data->output_fp) {
            pthread_mutex_lock(data->output_mutex);
            int write_status = fprintf(data->output_fp, ">%s:%"PRIhts_pos"-%"PRIhts_pos"\n%s\n", 
                   seq_name, start, end, seq);
            if (write_status < 0) {
                fprintf(stderr, "Thread %d: Error writing to output file: %s\n", 
                       data->thread_id, strerror(errno));
            }
            fflush(data->output_fp); // Ensure immediate write to disk
            pthread_mutex_unlock(data->output_mutex);
            
            if (data->config->verbose) {
                printf("Thread %d: Wrote sequence %s:%"PRIhts_pos"-%"PRIhts_pos" to output file\n",
                      data->thread_id, seq_name, start, end);
            }
        }
        
        free(seq);
    }
    
    // End timing
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    data->elapsed_time = (end_time.tv_sec - start_time.tv_sec) + 
                        (end_time.tv_nsec - start_time.tv_nsec) / 1e9;
    data->num_bases = bases_fetched;
    
    if (data->config->verbose) {
        printf("Thread %d: Fetched %"PRIu64" bases in %.3f seconds (%.2f bases/sec)\n", 
               data->thread_id, bases_fetched, data->elapsed_time, 
               bases_fetched / data->elapsed_time);
    }
    
    faidx_reader_destroy(reader);
    return NULL;
}

int main(int argc, char **argv) {
    bench_config_t config;
    faidx_meta_t *meta;
    pthread_t *threads;
    thread_data_t *thread_data;
    pthread_mutex_t output_mutex;
    FILE *output_fp = NULL;
    int i;
    double total_time = 0.0;
    uint64_t total_bases = 0;
    
    // Parse command line arguments
    config = parse_args(argc, argv);
    
    printf("Benchmark configuration:\n");
    printf("  FASTA file:  %s\n", config.fasta_file);
    printf("  Threads:     %d\n", config.num_threads);
    printf("  Seq count:   %d per thread (%d total)\n", config.seq_count, 
           config.seq_count * config.num_threads);
    printf("  Seq length:  %d\n", config.seq_length);
    printf("  Output:      %s\n", config.output_file ? config.output_file : "none");
    printf("  Seed:        %u\n", config.seed);
    printf("  Verbose:     %s\n", config.verbose ? "yes" : "no");
    
    // Output file handling message
    if (!config.output_file) {
        if (config.verbose) {
            printf("\nNote: No output file specified. Use -o option to write sequences to a file.\n");
        }
    } else {
        printf("Output file: %s will be created/overwritten\n", config.output_file);
    }
    
    // Load the FASTA index metadata
    meta = faidx_meta_load(config.fasta_file, FAI_FASTA, FAI_CREATE);
    if (!meta) {
        fprintf(stderr, "Failed to load FASTA index\n");
        return 1;
    }
    
    printf("Loaded index with %d sequences\n", faidx_meta_nseq(meta));
    
    // Open output file if specified
    if (config.output_file) {
        output_fp = fopen(config.output_file, "w");
        if (!output_fp) {
            fprintf(stderr, "Failed to open output file: %s\n", config.output_file);
            faidx_meta_destroy(meta);
            return 1;
        }
        pthread_mutex_init(&output_mutex, NULL);
    }
    
    // Allocate thread resources
    threads = malloc(config.num_threads * sizeof(pthread_t));
    thread_data = malloc(config.num_threads * sizeof(thread_data_t));
    if (!threads || !thread_data) {
        fprintf(stderr, "Failed to allocate memory for threads\n");
        if (output_fp) fclose(output_fp);
        faidx_meta_destroy(meta);
        return 1;
    }
    
    // Create and start threads
    for (i = 0; i < config.num_threads; i++) {
        thread_data[i].thread_id = i;
        thread_data[i].meta = meta;
        thread_data[i].config = &config;
        thread_data[i].seed = config.seed;
        thread_data[i].output_mutex = &output_mutex;
        thread_data[i].output_fp = output_fp;
        thread_data[i].num_bases = 0;
        thread_data[i].elapsed_time = 0.0;
        
        if (pthread_create(&threads[i], NULL, worker_thread, &thread_data[i]) != 0) {
            fprintf(stderr, "Failed to create thread %d\n", i);
            faidx_meta_destroy(meta);
            free(threads);
            free(thread_data);
            if (output_fp) {
                fclose(output_fp);
                pthread_mutex_destroy(&output_mutex);
            }
            return 1;
        }
    }
    
    // Wait for all threads to finish
    for (i = 0; i < config.num_threads; i++) {
        pthread_join(threads[i], NULL);
        total_time += thread_data[i].elapsed_time;
        total_bases += thread_data[i].num_bases;
    }
    
    // Calculate average time and throughput
    double avg_time = total_time / config.num_threads;
    double throughput = total_bases / avg_time;
    
    printf("\nBenchmark Results:\n");
    printf("  Total sequences fetched: %d\n", config.seq_count * config.num_threads);
    printf("  Total bases fetched:     %"PRIu64"\n", total_bases);
    printf("  Average time per thread: %.3f seconds\n", avg_time);
    printf("  Total throughput:        %.2f bases/second\n", throughput);
    
    // Clean up
    free(threads);
    free(thread_data);
    faidx_meta_destroy(meta);
    
    if (output_fp) {
        fclose(output_fp);
        pthread_mutex_destroy(&output_mutex);
        
        // Check if output file was created successfully
        FILE *check_fp = fopen(config.output_file, "r");
        if (check_fp) {
            fseek(check_fp, 0, SEEK_END);
            long file_size = ftell(check_fp);
            fclose(check_fp);
            printf("Sequences written to %s (size: %ld bytes)\n", config.output_file, file_size);
        } else {
            fprintf(stderr, "Warning: Cannot verify output file %s: %s\n", 
                   config.output_file, strerror(errno));
        }
    }
    
    return 0;
}
