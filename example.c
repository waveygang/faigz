#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <inttypes.h>

/* Include the implementation */
#define REENTRANT_FAIDX_IMPLEMENTATION
#include "faigz.h"

#define NUM_THREADS 4
#define NUM_FETCHES 10

typedef struct {
    faidx_meta_t *meta;
    int thread_id;
} thread_data_t;

void *worker_thread(void *arg) {
    thread_data_t *data = (thread_data_t *)arg;
    faidx_reader_t *reader;
    char *seq;
    hts_pos_t len;
    int i;
    
    printf("Thread %d: Creating reader\n", data->thread_id);
    reader = faidx_reader_create(data->meta);
    if (!reader) {
        fprintf(stderr, "Thread %d: Failed to create reader\n", data->thread_id);
        return NULL;
    }
    
    /* Get first sequence name */
    const char *first_seq = faidx_meta_iseq(data->meta, 0);
    if (!first_seq) {
        fprintf(stderr, "Thread %d: No sequences found\n", data->thread_id);
        faidx_reader_destroy(reader);
        return NULL;
    }
    
    /* Get sequence length */
    hts_pos_t seq_len = faidx_meta_seq_len(data->meta, first_seq);
    if (seq_len <= 0) {
        fprintf(stderr, "Thread %d: Invalid sequence length\n", data->thread_id);
        faidx_reader_destroy(reader);
        return NULL;
    }
    
    /* Fetch different parts of the sequence */
    for (i = 0; i < NUM_FETCHES; i++) {
        hts_pos_t start = (seq_len / NUM_FETCHES) * i;
        hts_pos_t end = start + 100;  /* Get 100 bases */
        
        if (end > seq_len) end = seq_len;
        
        seq = faidx_reader_fetch_seq(reader, first_seq, start, end - 1, &len);
        if (!seq) {
            fprintf(stderr, "Thread %d: Failed to fetch sequence\n", data->thread_id);
            continue;
        }
        
        printf("Thread %d: Fetch %d: %s:%"PRIhts_pos"-%"PRIhts_pos" length=%"PRIhts_pos" data=%.*s\n",
               data->thread_id, i, first_seq, start, end - 1, len, len > 20 ? 20 : (int)len, seq);
        
        free(seq);
    }
    
    printf("Thread %d: Cleaning up\n", data->thread_id);
    faidx_reader_destroy(reader);
    return NULL;
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <fasta_file>\n", argv[0]);
        return 1;
    }
    
    /* Load the index metadata once */
    faidx_meta_t *meta = faidx_meta_load(argv[1], FAI_FASTA, FAI_CREATE);
    if (!meta) {
        fprintf(stderr, "Failed to load FASTA index\n");
        return 1;
    }
    
    printf("Loaded index with %d sequences\n", faidx_meta_nseq(meta));
    
    /* Create worker threads */
    pthread_t threads[NUM_THREADS];
    thread_data_t thread_data[NUM_THREADS];
    
    for (int i = 0; i < NUM_THREADS; i++) {
        thread_data[i].meta = meta;
        thread_data[i].thread_id = i;
        
        if (pthread_create(&threads[i], NULL, worker_thread, &thread_data[i]) != 0) {
            fprintf(stderr, "Failed to create thread %d\n", i);
            faidx_meta_destroy(meta);
            return 1;
        }
    }
    
    /* Wait for threads to finish */
    for (int i = 0; i < NUM_THREADS; i++) {
        pthread_join(threads[i], NULL);
    }
    
    /* Clean up */
    faidx_meta_destroy(meta);
    
    printf("All done!\n");
    return 0;
}
