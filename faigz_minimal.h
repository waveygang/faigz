#ifndef FAIGZ_MINIMAL_H
#define FAIGZ_MINIMAL_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <pthread.h>
#include <inttypes.h>
#include <zlib.h>

#ifdef __cplusplus
extern "C" {
#endif

// Forward declarations
typedef struct faidx_meta_t faidx_meta_t;
typedef struct faidx_reader_t faidx_reader_t;

// Format options
typedef enum {
    FAI_NONE = 0,
    FAI_FASTA = 1,
    FAI_FASTQ = 2
} fai_format_options;

// Flags for faidx_meta_load
#define FAI_CREATE 0x01

// Position type
typedef int64_t hts_pos_t;

// Index entry structure
typedef struct {
    int id;
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} faidx1_t;

// Simple string hash table
typedef struct {
    char *key;
    faidx1_t val;
} hash_entry_t;

typedef struct {
    hash_entry_t *entries;
    int n_entries;
    int capacity;
} simple_hash_t;

// Shared metadata structure
struct faidx_meta_t {
    int n, m;                     // Sequence count and allocation size
    char **name;                  // Array of sequence names
    simple_hash_t *hash;          // Hash table mapping names to positions
    fai_format_options format;    // FAI_FASTA or FAI_FASTQ
    
    // Source file paths
    char *fasta_path;            // Path to the FASTA/FASTQ file
    char *fai_path;              // Path to the .fai index
    char *gzi_path;              // Path to the .gzi index (if using BGZF)
    
    // Reference count and mutex for thread safety
    int ref_count;
    pthread_mutex_t mutex;
    
    // Flag indicating if the source is BGZF compressed
    int is_bgzf;
};

// Reader structure containing thread-specific data
struct faidx_reader_t {
    faidx_meta_t *meta;          // Shared metadata (not owned)
    FILE *fp;                    // File pointer for reading
    gzFile gzfp;                 // gzFile pointer for compressed files
};

// Function declarations
faidx_meta_t *faidx_meta_load(const char *filename, fai_format_options format, int flags);
faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta);
void faidx_meta_destroy(faidx_meta_t *meta);
faidx_reader_t *faidx_reader_create(faidx_meta_t *meta);
void faidx_reader_destroy(faidx_reader_t *reader);
char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name,
                           hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);
char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name,
                            hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);
int faidx_meta_nseq(const faidx_meta_t *meta);
const char *faidx_meta_iseq(const faidx_meta_t *meta, int i);
hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq);
int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq);

#ifdef __cplusplus
}
#endif

#endif /* FAIGZ_MINIMAL_H */