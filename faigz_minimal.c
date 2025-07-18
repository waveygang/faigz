#include "faigz_minimal.h"
#include <ctype.h>
#include <sys/stat.h>
#include <unistd.h>

// Hash table implementation
static simple_hash_t *hash_init(void) {
    simple_hash_t *h = calloc(1, sizeof(simple_hash_t));
    if (!h) return NULL;
    h->capacity = 16;
    h->entries = calloc(h->capacity, sizeof(hash_entry_t));
    if (!h->entries) {
        free(h);
        return NULL;
    }
    return h;
}

static void hash_destroy(simple_hash_t *h) {
    if (!h) return;
    for (int i = 0; i < h->n_entries; i++) {
        free(h->entries[i].key);
    }
    free(h->entries);
    free(h);
}

static int hash_put(simple_hash_t *h, const char *key, faidx1_t val) {
    if (h->n_entries >= h->capacity) {
        // Simple resize
        h->capacity *= 2;
        h->entries = realloc(h->entries, h->capacity * sizeof(hash_entry_t));
        if (!h->entries) return -1;
    }
    
    h->entries[h->n_entries].key = strdup(key);
    h->entries[h->n_entries].val = val;
    h->n_entries++;
    return 0;
}

static faidx1_t *hash_get(simple_hash_t *h, const char *key) {
    for (int i = 0; i < h->n_entries; i++) {
        if (strcmp(h->entries[i].key, key) == 0) {
            return &h->entries[i].val;
        }
    }
    return NULL;
}

// Utility functions
static char *str_dup(const char *str) {
    if (!str) return NULL;
    int len = strlen(str) + 1;
    char *s = malloc(len);
    if (!s) return NULL;
    memcpy(s, str, len);
    return s;
}

static int is_bgzf_file(const char *filename) {
    FILE *fp = fopen(filename, "rb");
    if (!fp) return 0;
    
    unsigned char magic[2];
    int is_bgzf = 0;
    if (fread(magic, 1, 2, fp) == 2) {
        is_bgzf = (magic[0] == 0x1f && magic[1] == 0x8b);
    }
    fclose(fp);
    return is_bgzf;
}

static int create_fai_index(const char *fasta_path, const char *fai_path) {
    FILE *fasta_fp = fopen(fasta_path, "r");
    if (!fasta_fp) return -1;
    
    FILE *fai_fp = fopen(fai_path, "w");
    if (!fai_fp) {
        fclose(fasta_fp);
        return -1;
    }
    
    char line[1024];
    char seq_name[256] = {0};
    uint64_t seq_len = 0;
    uint64_t seq_offset = 0;
    uint32_t line_blen = 0;
    uint32_t line_len = 0;
    int in_sequence = 0;
    uint64_t current_offset = 0;
    
    while (fgets(line, sizeof(line), fasta_fp)) {
        int line_length = strlen(line);
        
        if (line[0] == '>') {
            // If we were processing a sequence, write its index entry
            if (in_sequence && seq_name[0]) {
                fprintf(fai_fp, "%s\t%lu\t%lu\t%u\t%u\n", 
                       seq_name, seq_len, seq_offset, line_blen, line_len);
            }
            
            // Start new sequence
            in_sequence = 1;
            seq_len = 0;
            seq_offset = current_offset + line_length;
            line_blen = 0;
            line_len = 0;
            
            // Extract sequence name (everything after '>' until whitespace)
            char *name_start = line + 1;
            char *name_end = name_start;
            while (*name_end && *name_end != ' ' && *name_end != '\t' && *name_end != '\n') {
                name_end++;
            }
            int name_len = name_end - name_start;
            if (name_len >= sizeof(seq_name)) name_len = sizeof(seq_name) - 1;
            strncpy(seq_name, name_start, name_len);
            seq_name[name_len] = '\0';
            
        } else if (in_sequence && line[0] != '\n' && line[0] != '\r') {
            // Count sequence characters (excluding newlines)
            int bases_in_line = 0;
            for (int i = 0; i < line_length; i++) {
                if (line[i] != '\n' && line[i] != '\r') {
                    bases_in_line++;
                }
            }
            
            seq_len += bases_in_line;
            
            // Set line length info from first line of sequence
            if (line_blen == 0) {
                line_blen = bases_in_line;
                line_len = line_length;
            }
        }
        
        current_offset += line_length;
    }
    
    // Write final sequence entry
    if (in_sequence && seq_name[0]) {
        fprintf(fai_fp, "%s\t%lu\t%lu\t%u\t%u\n", 
               seq_name, seq_len, seq_offset, line_blen, line_len);
    }
    
    fclose(fasta_fp);
    fclose(fai_fp);
    return 0;
}

static int load_fai_index(faidx_meta_t *meta, const char *fai_path) {
    FILE *fp = fopen(fai_path, "r");
    if (!fp) return -1;
    
    char line[1024];
    int idx = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        char *name = strtok(line, "\t");
        char *len_str = strtok(NULL, "\t");
        char *offset_str = strtok(NULL, "\t");
        char *line_blen_str = strtok(NULL, "\t");
        char *line_len_str = strtok(NULL, "\t\n");
        
        if (!name || !len_str || !offset_str || !line_blen_str || !line_len_str) {
            continue;
        }
        
        // Expand arrays if needed
        if (idx >= meta->m) {
            meta->m = meta->m ? meta->m * 2 : 16;
            meta->name = realloc(meta->name, meta->m * sizeof(char*));
            if (!meta->name) {
                fclose(fp);
                return -1;
            }
        }
        
        meta->name[idx] = str_dup(name);
        if (!meta->name[idx]) {
            fclose(fp);
            return -1;
        }
        
        faidx1_t val;
        val.id = idx;
        val.len = atoll(len_str);
        val.seq_offset = atoll(offset_str);
        val.line_blen = atoi(line_blen_str);
        val.line_len = atoi(line_len_str);
        val.qual_offset = 0; // For FASTQ, this would be calculated
        
        if (hash_put(meta->hash, name, val) < 0) {
            fclose(fp);
            return -1;
        }
        
        idx++;
    }
    
    meta->n = idx;
    fclose(fp);
    return 0;
}

// Public API implementation
faidx_meta_t *faidx_meta_load(const char *filename, fai_format_options format, int flags) {
    if (!filename) return NULL;
    
    faidx_meta_t *meta = calloc(1, sizeof(faidx_meta_t));
    if (!meta) return NULL;
    
    // Initialize mutex
    if (pthread_mutex_init(&meta->mutex, NULL) != 0) {
        free(meta);
        return NULL;
    }
    
    meta->format = format;
    meta->ref_count = 1;
    meta->is_bgzf = is_bgzf_file(filename);
    
    // Store file paths
    meta->fasta_path = str_dup(filename);
    
    // Construct index paths
    char fai_path[1024];
    snprintf(fai_path, sizeof(fai_path), "%s.fai", filename);
    meta->fai_path = str_dup(fai_path);
    
    char gzi_path[1024];
    snprintf(gzi_path, sizeof(gzi_path), "%s.gzi", filename);
    meta->gzi_path = str_dup(gzi_path);
    
    if (!meta->fasta_path || !meta->fai_path || !meta->gzi_path) {
        faidx_meta_destroy(meta);
        return NULL;
    }
    
    // Initialize hash table
    meta->hash = hash_init();
    if (!meta->hash) {
        faidx_meta_destroy(meta);
        return NULL;
    }
    
    // Try to load the index, or create it if it doesn't exist and FAI_CREATE is set
    if (load_fai_index(meta, meta->fai_path) < 0) {
        if (flags & FAI_CREATE) {
            // Try to create the index
            if (create_fai_index(meta->fasta_path, meta->fai_path) < 0) {
                faidx_meta_destroy(meta);
                return NULL;
            }
            // Now try to load the newly created index
            if (load_fai_index(meta, meta->fai_path) < 0) {
                faidx_meta_destroy(meta);
                return NULL;
            }
        } else {
            faidx_meta_destroy(meta);
            return NULL;
        }
    }
    
    return meta;
}

faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta) {
    if (!meta) return NULL;
    
    pthread_mutex_lock(&meta->mutex);
    meta->ref_count++;
    pthread_mutex_unlock(&meta->mutex);
    
    return meta;
}

void faidx_meta_destroy(faidx_meta_t *meta) {
    if (!meta) return;
    
    pthread_mutex_lock(&meta->mutex);
    meta->ref_count--;
    int should_free = (meta->ref_count <= 0);
    pthread_mutex_unlock(&meta->mutex);
    
    if (should_free) {
        if (meta->hash) hash_destroy(meta->hash);
        
        if (meta->name) {
            for (int i = 0; i < meta->n; i++) {
                free(meta->name[i]);
            }
            free(meta->name);
        }
        
        free(meta->fasta_path);
        free(meta->fai_path);
        free(meta->gzi_path);
        
        pthread_mutex_destroy(&meta->mutex);
        free(meta);
    }
}

faidx_reader_t *faidx_reader_create(faidx_meta_t *meta) {
    if (!meta) return NULL;
    
    faidx_reader_t *reader = calloc(1, sizeof(faidx_reader_t));
    if (!reader) return NULL;
    
    reader->meta = faidx_meta_ref(meta);
    
    // Open file
    if (meta->is_bgzf) {
        reader->gzfp = gzopen(meta->fasta_path, "r");
        if (!reader->gzfp) {
            faidx_meta_destroy(reader->meta);
            free(reader);
            return NULL;
        }
    } else {
        reader->fp = fopen(meta->fasta_path, "r");
        if (!reader->fp) {
            faidx_meta_destroy(reader->meta);
            free(reader);
            return NULL;
        }
    }
    
    return reader;
}

void faidx_reader_destroy(faidx_reader_t *reader) {
    if (!reader) return;
    
    if (reader->fp) fclose(reader->fp);
    if (reader->gzfp) gzclose(reader->gzfp);
    
    faidx_meta_destroy(reader->meta);
    free(reader);
}

char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name,
                           hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len) {
    if (!reader || !c_name) return NULL;
    
    faidx1_t *entry = hash_get(reader->meta->hash, c_name);
    if (!entry) return NULL;
    
    // Adjust coordinates
    if (p_beg_i < 0) p_beg_i = 0;
    if (p_end_i < 0 || p_end_i > entry->len) p_end_i = entry->len;
    if (p_beg_i >= p_end_i) return NULL;
    
    hts_pos_t seq_len = p_end_i - p_beg_i;
    char *seq = malloc(seq_len + 1);
    if (!seq) return NULL;
    
    // Simple implementation - seek to position and read
    // This is a simplified version that doesn't handle all edge cases
    if (reader->meta->is_bgzf) {
        gzseek(reader->gzfp, entry->seq_offset, SEEK_SET);
        // Skip to the right line and position
        // This is a simplified implementation
        if (len) *len = 0;
        free(seq);
        return NULL; // Not implemented for compressed files in this minimal version
    } else {
        fseek(reader->fp, entry->seq_offset, SEEK_SET);
        
        // Skip to the correct position
        hts_pos_t pos = 0;
        hts_pos_t line_pos = 0;
        int c;
        
        while (pos < p_beg_i && (c = fgetc(reader->fp)) != EOF) {
            if (c == '\n') {
                line_pos = 0;
            } else if (c != '\r') {
                pos++;
                line_pos++;
            }
        }
        
        // Read the sequence
        hts_pos_t read_len = 0;
        while (read_len < seq_len && (c = fgetc(reader->fp)) != EOF) {
            if (c == '\n' || c == '\r') {
                continue;
            }
            if (c == '>' || c == '+') {
                break; // Hit next sequence
            }
            seq[read_len++] = c;
        }
        
        seq[read_len] = '\0';
        if (len) *len = read_len;
        return seq;
    }
}

char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name,
                            hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len) {
    if (!reader || reader->meta->format != FAI_FASTQ) return NULL;
    // Quality string fetching not implemented in this minimal version
    if (len) *len = 0;
    return NULL;
}

int faidx_meta_nseq(const faidx_meta_t *meta) {
    return meta ? meta->n : 0;
}

const char *faidx_meta_iseq(const faidx_meta_t *meta, int i) {
    return (meta && i >= 0 && i < meta->n) ? meta->name[i] : NULL;
}

hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq) {
    if (!meta || !seq) return -1;
    
    faidx1_t *entry = hash_get(meta->hash, seq);
    return entry ? entry->len : -1;
}

int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq) {
    if (!meta || !seq) return 0;
    
    faidx1_t *entry = hash_get(meta->hash, seq);
    return entry != NULL;
}