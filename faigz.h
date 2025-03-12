#ifndef REENTRANT_FAIDX_H
#define REENTRANT_FAIDX_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <pthread.h>

#include "htslib/bgzf.h"
#include "htslib/faidx.h"
#include "htslib/khash.h"
#include "htslib/hts.h"

#ifdef __cplusplus
extern "C" {
#endif

// Key structures needed for our implementation
typedef struct faidx1_t {
    int id;
    uint32_t line_len, line_blen;
    uint64_t len;
    uint64_t seq_offset;
    uint64_t qual_offset;
} faidx1_t;

KHASH_MAP_INIT_STR(s, faidx1_t)

// Shared metadata structure containing only the indices
typedef struct {
    int n, m;                     // Sequence count and allocation size
    char **name;                  // Array of sequence names
    khash_t(s) *hash;            // Hash table mapping names to positions
    enum fai_format_options format; // FAI_FASTA or FAI_FASTQ
    
    // Source file paths
    char *fasta_path;            // Path to the FASTA/FASTQ file
    char *fai_path;              // Path to the .fai index
    char *gzi_path;              // Path to the .gzi index (if using BGZF)
    
    // BGZF index (if applicable)
    struct bgzidx_t *bgzf_idx;   // Shared BGZF index
    
    // Reference count and mutex for thread safety
    int ref_count;
    pthread_mutex_t mutex;
    
    // Flag indicating if the source is BGZF compressed
    int is_bgzf;
} faidx_meta_t;

// Reader structure containing thread-specific data
typedef struct {
    faidx_meta_t *meta;          // Shared metadata (not owned)
    BGZF *bgzf;                  // Thread-local BGZF handle
} faidx_reader_t;

/**
 * Load FASTA/FASTQ index metadata without opening the file.
 * 
 * @param filename Path to the FASTA/FASTQ file
 * @param format FAI_FASTA or FAI_FASTQ
 * @param flags Option flags (see FAI_CREATE in faidx.h)
 * @return Pointer to metadata or NULL on error
 */
faidx_meta_t *faidx_meta_load(const char *filename, enum fai_format_options format, int flags);

/**
 * Increment reference count on metadata
 * 
 * @param meta Metadata to reference
 * @return The metadata
 */
faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta);

/**
 * Free shared metadata and its resources.
 * Decrements the reference count and frees only when it reaches zero.
 * 
 * @param meta Metadata to free
 */
void faidx_meta_destroy(faidx_meta_t *meta);

/**
 * Create a reader from shared metadata
 * 
 * @param meta Shared metadata (reference count is incremented)
 * @return New reader or NULL on error
 */
faidx_reader_t *faidx_reader_create(faidx_meta_t *meta);

/**
 * Destroy a reader.
 * This does not affect the shared metadata.
 * 
 * @param reader Reader to destroy
 */
void faidx_reader_destroy(faidx_reader_t *reader);

/**
 * Fetch sequence from a specific region
 * 
 * @param reader Reader to use
 * @param c_name Region name
 * @param p_beg_i Beginning position (0-based)
 * @param p_end_i End position (0-based)
 * @param len Output parameter for sequence length
 * @return Sequence string (must be freed by caller) or NULL on error
 */
char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name,
                           hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);

/**
 * Fetch the quality string for a specific region (FASTQ only)
 * 
 * @param reader Reader to use
 * @param c_name Region name
 * @param p_beg_i Beginning position (0-based)
 * @param p_end_i End position (0-based)
 * @param len Output parameter for string length
 * @return Quality string (must be freed by caller) or NULL on error
 */
char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name,
                            hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len);

/**
 * Get number of sequences in the index
 * 
 * @param meta Metadata
 * @return Number of sequences
 */
int faidx_meta_nseq(const faidx_meta_t *meta);

/**
 * Get name of the i-th sequence
 * 
 * @param meta Metadata
 * @param i Sequence index
 * @return Sequence name or NULL
 */
const char *faidx_meta_iseq(const faidx_meta_t *meta, int i);

/**
 * Get sequence length
 * 
 * @param meta Metadata
 * @param seq Sequence name
 * @return Sequence length or -1 if not found
 */
hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq);

/**
 * Check if a sequence exists in the index
 * 
 * @param meta Metadata
 * @param seq Sequence name
 * @return 1 if present, 0 if absent
 */
int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq);

/**
 * Parse a region string
 * 
 * @param meta Metadata
 * @param s Region string
 * @param tid Output parameter for sequence ID
 * @param beg Output parameter for beginning position
 * @param end Output parameter for end position
 * @param flags Parsing flags
 * @return Pointer to end of parsed region or NULL on error
 */
const char *faidx_meta_parse_region(const faidx_meta_t *meta, const char *s,
                                  int *tid, hts_pos_t *beg, hts_pos_t *end,
                                  int flags);

#ifdef __cplusplus
}
#endif

/* Internal declarations for implementation */
#ifdef REENTRANT_FAIDX_IMPLEMENTATION

// Helper functions for internal use
static char *kstrdup(const char *str);
static faidx1_t *fai_get_val(const faidx_meta_t *meta, const char *str,
                           hts_pos_t *len, faidx1_t *val, hts_pos_t *fbeg, hts_pos_t *fend);
static int faidx_adjust_position(const faidx_meta_t *meta, int end_adjust,
                              faidx1_t *val_out, const char *c_name,
                              hts_pos_t *p_beg_i, hts_pos_t *p_end_i,
                              hts_pos_t *len);
static int fai_name2id(void *v, const char *ref);
static char *fai_retrieve(BGZF *bgzf, const faidx1_t *val,
                       uint64_t offset, hts_pos_t beg, hts_pos_t end, hts_pos_t *len);
static struct bgzidx_t *bgzf_index_load_direct(const char *bname);
static BGZF *bgzf_open_shared_idx(const char *path, const char *mode, struct bgzidx_t *idx);

/* Implementation of the functions */

static char *kstrdup(const char *str) {
    int len = strlen(str) + 1;
    char *s = (char*)malloc(len);
    if (!s) return NULL;
    memcpy(s, str, len);
    return s;
}

static int fai_name2id(void *v, const char *ref) {
    faidx_meta_t *meta = (faidx_meta_t*)v;
    khint_t k = kh_get(s, meta->hash, ref);
    return k == kh_end(meta->hash) ? -1 : kh_val(meta->hash, k).id;
}

/* This function loads the BGZF index directly without opening the file */
static struct bgzidx_t *bgzf_index_load_direct(const char *bname) {
    struct bgzidx_t *idx = NULL;
    hFILE *fp = NULL;
    
    fp = hopen(bname, "rb");
    if (!fp) return NULL;
    
    idx = (struct bgzidx_t*)calloc(1, sizeof(struct bgzidx_t));
    if (!idx) {
        hclose(fp);
        return NULL;
    }
    
    uint64_t x;
    if (hread_uint64(&x, fp) < 0) goto fail;
    
    idx->noffs = idx->moffs = x + 1;
    idx->offs = (bgzidx1_t*)malloc(idx->moffs * sizeof(bgzidx1_t));
    if (!idx->offs) goto fail;
    idx->offs[0].caddr = idx->offs[0].uaddr = 0;
    
    int i;
    for (i = 1; i < idx->noffs; i++) {
        if (hread_uint64(&idx->offs[i].caddr, fp) < 0) goto fail;
        if (hread_uint64(&idx->offs[i].uaddr, fp) < 0) goto fail;
    }
    
    hclose(fp);
    return idx;
    
fail:
    if (idx) {
        free(idx->offs);
        free(idx);
    }
    if (fp) hclose(fp);
    return NULL;
}

/* This function opens a BGZF file with a pre-loaded index */
static BGZF *bgzf_open_shared_idx(const char *path, const char *mode, struct bgzidx_t *idx) {
    BGZF *bgzf = bgzf_open(path, mode);
    if (!bgzf) return NULL;
    
    /* Clean up any index that might have been loaded */
    if (bgzf->idx) {
        free(bgzf->idx->offs);
        free(bgzf->idx);
    }
    
    /* Set the shared index */
    bgzf->idx = idx;
    
    return bgzf;
}

/* Load metadata from a FASTA/FASTQ file */
faidx_meta_t *faidx_meta_load(const char *filename, enum fai_format_options format, int flags) {
    kstring_t fai_kstr = {0}, gzi_kstr = {0};
    faidx_t *fai = NULL;
    faidx_meta_t *meta = NULL;
    
    /* Handle NULL filename */
    if (!filename) {
        errno = EINVAL;
        return NULL;
    }
    
    /* Construct file paths */
    const char *fnfai = NULL, *fngzi = NULL;
    
    if (ksprintf(&fai_kstr, "%s.fai", filename) < 0) goto fail;
    fnfai = fai_kstr.s;
    
    if (ksprintf(&gzi_kstr, "%s.gzi", filename) < 0) goto fail;
    fngzi = gzi_kstr.s;
    
    /* Load the FASTA/FASTQ index */
    fai = fai_load3_format(filename, fnfai, fngzi, flags, format);
    if (!fai) goto fail;
    
    /* Create the metadata structure */
    meta = (faidx_meta_t*)calloc(1, sizeof(faidx_meta_t));
    if (!meta) goto fail;
    
    /* Initialize the mutex */
    if (pthread_mutex_init(&meta->mutex, NULL) != 0) goto fail;
    
    /* Copy the FAI data */
    meta->n = fai->n;
    meta->m = fai->m;
    meta->format = fai->format;
    meta->ref_count = 1;
    
    /* Copy sequence names */
    meta->name = (char**)malloc(meta->m * sizeof(char*));
    if (!meta->name) goto fail;
    
    for (int i = 0; i < meta->n; i++) {
        meta->name[i] = kstrdup(fai->name[i]);
        if (!meta->name[i]) goto fail;
    }
    
    /* Create and populate hash table */
    meta->hash = kh_init(s);
    if (!meta->hash) goto fail;
    
    for (int i = 0; i < meta->n; i++) {
        int absent;
        khint_t k = kh_put(s, meta->hash, meta->name[i], &absent);
        if (absent) {
            khint_t fai_k = kh_get(s, fai->hash, fai->name[i]);
            kh_val(meta->hash, k) = kh_val(fai->hash, fai_k);
        }
    }
    
    /* Store file paths */
    meta->fasta_path = kstrdup(filename);
    meta->fai_path = kstrdup(fnfai);
    meta->gzi_path = kstrdup(fngzi);
    if (!meta->fasta_path || !meta->fai_path || !meta->gzi_path) goto fail;
    
    /* Check if file is BGZF compressed and load BGZF index if needed */
    meta->is_bgzf = fai->bgzf && fai->bgzf->is_compressed;
    
    if (meta->is_bgzf) {
        /* We need to load the BGZF index separately */
        hFILE *fp = hopen(fngzi, "rb");
        if (fp) {
            hclose(fp); /* Just checking existence */
            meta->bgzf_idx = bgzf_index_load_direct(fngzi);
        }
    }
    
    /* Clean up */
    free(fai_kstr.s);
    free(gzi_kstr.s);
    fai_destroy(fai);
    
    return meta;
    
fail:
    if (meta) {
        if (meta->hash) kh_destroy(s, meta->hash);
        if (meta->name) {
            for (int i = 0; i < meta->n; i++) free(meta->name[i]);
            free(meta->name);
        }
        free(meta->fasta_path);
        free(meta->fai_path);
        free(meta->gzi_path);
        if (meta->bgzf_idx) {
            free(meta->bgzf_idx->offs);
            free(meta->bgzf_idx);
        }
        pthread_mutex_destroy(&meta->mutex);
        free(meta);
    }
    free(fai_kstr.s);
    free(gzi_kstr.s);
    if (fai) fai_destroy(fai);
    return NULL;
}

/* Reference counting for metadata */
faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta) {
    if (!meta) return NULL;
    
    pthread_mutex_lock(&meta->mutex);
    meta->ref_count++;
    pthread_mutex_unlock(&meta->mutex);
    
    return meta;
}

/* Destroy metadata */
void faidx_meta_destroy(faidx_meta_t *meta) {
    if (!meta) return;
    
    pthread_mutex_lock(&meta->mutex);
    meta->ref_count--;
    int should_free = (meta->ref_count <= 0);
    pthread_mutex_unlock(&meta->mutex);
    
    if (should_free) {
        /* Free all resources */
        if (meta->hash) {
            kh_destroy(s, meta->hash);
        }
        
        if (meta->name) {
            for (int i = 0; i < meta->n; i++) {
                free(meta->name[i]);
            }
            free(meta->name);
        }
        
        if (meta->bgzf_idx) {
            free(meta->bgzf_idx->offs);
            free(meta->bgzf_idx);
        }
        
        free(meta->fasta_path);
        free(meta->fai_path);
        free(meta->gzi_path);
        
        pthread_mutex_destroy(&meta->mutex);
        free(meta);
    }
}

/* Create a reader */
faidx_reader_t *faidx_reader_create(faidx_meta_t *meta) {
    if (!meta) return NULL;
    
    faidx_reader_t *reader = (faidx_reader_t*)calloc(1, sizeof(faidx_reader_t));
    if (!reader) return NULL;
    
    /* Reference the metadata */
    reader->meta = faidx_meta_ref(meta);
    
    /* Open the BGZF file with shared index */
    if (meta->is_bgzf) {
        reader->bgzf = bgzf_open_shared_idx(meta->fasta_path, "rb", meta->bgzf_idx);
    } else {
        reader->bgzf = bgzf_open(meta->fasta_path, "rb");
    }
    
    if (!reader->bgzf) {
        faidx_meta_destroy(reader->meta);
        free(reader);
        return NULL;
    }
    
    return reader;
}

/* Destroy a reader */
void faidx_reader_destroy(faidx_reader_t *reader) {
    if (!reader) return;
    
    if (reader->bgzf) {
        /* Don't destroy the shared BGZF index */
        reader->bgzf->idx = NULL;
        bgzf_close(reader->bgzf);
    }
    
    faidx_meta_destroy(reader->meta);
    free(reader);
}

/* Helper: Get the faidx1_t value for a region */
static faidx1_t *fai_get_val(const faidx_meta_t *meta, const char *str,
                           hts_pos_t *len, faidx1_t *val, hts_pos_t *fbeg, hts_pos_t *fend) {
    khiter_t iter;
    int id;
    hts_pos_t beg, end;
    
    if (!faidx_meta_parse_region(meta, str, &id, &beg, &end, 0)) {
        if (len) *len = -2;
        return NULL;
    }
    
    khash_t(s) *h = meta->hash;
    iter = kh_get(s, h, meta->name[id]);
    if (iter >= kh_end(h)) {
        if (len) *len = -2;
        return NULL;
    }
    
    *val = kh_val(h, iter);
    
    if (beg >= val->len) beg = val->len;
    if (end >= val->len) end = val->len;
    if (beg > end) beg = end;
    
    *fbeg = beg;
    *fend = end;
    
    return val;
}

/* Helper: Adjust position to sequence boundaries */
static int faidx_adjust_position(const faidx_meta_t *meta, int end_adjust,
                              faidx1_t *val_out, const char *c_name,
                              hts_pos_t *p_beg_i, hts_pos_t *p_end_i,
                              hts_pos_t *len) {
    khiter_t iter;
    faidx1_t *val;
    
    /* Adjust position */
    iter = kh_get(s, meta->hash, c_name);
    
    if (iter == kh_end(meta->hash)) {
        if (len) *len = -2;
        return 1;
    }
    
    val = &kh_val(meta->hash, iter);
    
    if (val_out) *val_out = *val;
    
    if (*p_end_i < *p_beg_i) *p_beg_i = *p_end_i;
    
    if (*p_beg_i < 0) *p_beg_i = 0;
    else if (val->len <= *p_beg_i) *p_beg_i = val->len;
    
    if (*p_end_i < 0) *p_end_i = 0;
    else if (val->len <= *p_end_i) *p_end_i = val->len - end_adjust;
    
    return 0;
}

/* Helper: Retrieve sequence data */
static char *fai_retrieve(BGZF *bgzf, const faidx1_t *val,
                       uint64_t offset, hts_pos_t beg, hts_pos_t end, hts_pos_t *len) {
    char *buffer, *s;
    ssize_t nread, remaining, firstline_len, firstline_blen;
    int ret;
    
    if ((uint64_t)end - (uint64_t)beg >= SIZE_MAX - 2) {
        if (len) *len = -1;
        return NULL;
    }
    
    if (val->line_blen <= 0) {
        if (len) *len = -1;
        return NULL;
    }
    
    ret = bgzf_useek(bgzf, offset + beg / val->line_blen * val->line_len + beg % val->line_blen, SEEK_SET);
    
    if (ret < 0) {
        if (len) *len = -1;
        return NULL;
    }
    
    /* Over-allocate for one end-of-line sequence */
    buffer = (char*)malloc((size_t)end - beg + val->line_len - val->line_blen + 1);
    if (!buffer) {
        if (len) *len = -1;
        return NULL;
    }
    
    remaining = *len = end - beg;
    firstline_blen = val->line_blen - beg % val->line_blen;
    
    /* Special case when entire interval is within a single line */
    if (remaining <= firstline_blen) {
        nread = bgzf_read(bgzf, buffer, remaining);
        if (nread < remaining) {
            free(buffer);
            if (len) *len = -1;
            return NULL;
        }
        buffer[nread] = '\0';
        return buffer;
    }
    
    s = buffer;
    firstline_len = val->line_len - beg % val->line_blen;
    
    /* Read first line */
    nread = bgzf_read(bgzf, s, firstline_len);
    if (nread < firstline_len) {
        free(buffer);
        if (len) *len = -1;
        return NULL;
    }
    s += firstline_blen;
    remaining -= firstline_blen;
    
    /* Read complete lines */
    while (remaining > val->line_blen) {
        nread = bgzf_read(bgzf, s, val->line_len);
        if (nread < (ssize_t)val->line_len) {
            free(buffer);
            if (len) *len = -1;
            return NULL;
        }
        s += val->line_blen;
        remaining -= val->line_blen;
    }
    
    /* Read final partial line */
    if (remaining > 0) {
        nread = bgzf_read(bgzf, s, remaining);
        if (nread < remaining) {
            free(buffer);
            if (len) *len = -1;
            return NULL;
        }
        s += remaining;
    }
    
    *s = '\0';
    return buffer;
}

/* Fetch sequence */
char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name,
                          hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len) {
    faidx1_t val;
    
    /* Adjust position */
    if (faidx_adjust_position(reader->meta, 1, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }
    
    /* Retrieve sequence */
    return fai_retrieve(reader->bgzf, &val, val.seq_offset, p_beg_i, p_end_i + 1, len);
}

/* Fetch quality string */
char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name,
                            hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len) {
    faidx1_t val;
    
    if (reader->meta->format != FAI_FASTQ) {
        if (len) *len = -2;
        return NULL;
    }
    
    /* Adjust position */
    if (faidx_adjust_position(reader->meta, 1, &val, c_name, &p_beg_i, &p_end_i, len)) {
        return NULL;
    }
    
    /* Retrieve quality string */
    return fai_retrieve(reader->bgzf, &val, val.qual_offset, p_beg_i, p_end_i + 1, len);
}

/* Get number of sequences */
int faidx_meta_nseq(const faidx_meta_t *meta) {
    return meta ? meta->n : 0;
}

/* Get sequence name */
const char *faidx_meta_iseq(const faidx_meta_t *meta, int i) {
    return (meta && i >= 0 && i < meta->n) ? meta->name[i] : NULL;
}

/* Get sequence length */
hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq) {
    if (!meta || !seq) return -1;
    
    khint_t k = kh_get(s, meta->hash, seq);
    if (k == kh_end(meta->hash)) return -1;
    
    return kh_val(meta->hash, k).len;
}

/* Check if sequence exists */
int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq) {
    if (!meta || !seq) return 0;
    
    khint_t k = kh_get(s, meta->hash, seq);
    return (k != kh_end(meta->hash));
}

/* Parse a region string */
const char *faidx_meta_parse_region(const faidx_meta_t *meta, const char *s,
                                 int *tid, hts_pos_t *beg, hts_pos_t *end,
                                 int flags) {
    return hts_parse_region(s, tid, beg, end, fai_name2id, (void*)meta, flags);
}

#endif /* REENTRANT_FAIDX_IMPLEMENTATION */

#endif /* REENTRANT_FAIDX_H */
