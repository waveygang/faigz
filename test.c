#define REENTRANT_FAIDX_IMPLEMENTATION
#include "faigz.h"
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <assert.h>
#include <inttypes.h>

#define NUM_THREADS 4

typedef struct {
    faidx_meta_t *meta;
    int thread_id;
    int num_seqs;
    int success;
} thread_data_t;

void *test_thread(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    faidx_reader_t *reader;
    
    reader = faidx_reader_create(data->meta);
    if (!reader) {
        fprintf(stderr, "Thread %d: Failed to create reader\n", data->thread_id);
        return NULL;
    }
    
    // Test basic functionality for each sequence
    for (int i = 0; i < data->num_seqs; i++) {
        const char *seq_name = faidx_meta_iseq(data->meta, i);
        hts_pos_t seq_len = faidx_meta_seq_len(data->meta, seq_name);
        
        if (seq_len <= 0) continue;
        
        // Fetch first 10 bases (or less if shorter sequence)
        hts_pos_t fetch_len = seq_len < 10 ? seq_len : 10;
        hts_pos_t result_len;
        
        char *seq = faidx_reader_fetch_seq(reader, seq_name, 0, fetch_len - 1, &result_len);
        
        if (!seq) {
            fprintf(stderr, "Thread %d: Failed to fetch %s\n", data->thread_id, seq_name);
            continue;
        }
        
        if (result_len != fetch_len) {
            fprintf(stderr, "Thread %d: Length mismatch for %s: %"PRIhts_pos" vs %"PRIhts_pos"\n", 
                   data->thread_id, seq_name, result_len, fetch_len);
            free(seq);
            continue;
        }
        
        free(seq);
    }
    
    data->success = 1;
    faidx_reader_destroy(reader);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fasta_file>\n", argv[0]);
        return 1;
    }
    
    printf("Testing faigz library...\n");
    
    // Load index metadata
    faidx_meta_t *meta = faidx_meta_load(argv[1], FAI_FASTA, FAI_CREATE);
    if (!meta) {
        fprintf(stderr, "Failed to load FASTA index\n");
        return 1;
    }
    
    int num_seqs = faidx_meta_nseq(meta);
    printf("Loaded index with %d sequences\n", num_seqs);
    
    // Create worker threads
    pthread_t threads[NUM_THREADS];
    thread_data_t thread_data[NUM_THREADS];
    
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].meta = meta;
        thread_data[i].thread_id = i;
        thread_data[i].num_seqs = num_seqs;
        thread_data[i].success = 0;
        
        if (pthread_create(&threads[i], NULL, test_thread, &thread_data[i]) != 0) {
            fprintf(stderr, "Failed to create thread %d\n", i);
            faidx_meta_destroy(meta);
            return 1;
        }
    }
    
    // Wait for threads to finish
    int all_success = 1;
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
        all_success &= thread_data[i].success;
    }
    
    // Clean up
    faidx_meta_destroy(meta);
    
    if (all_success) {
        printf("All tests passed!\n");
        return 0;
    } else {
        fprintf(stderr, "Some tests failed\n");
        return 1;
    }
}
