# faigz - Reentrant FASTA/BGZF Index Library

A reentrant implementation of FASTA/BGZF indexing that efficiently shares both FASTA and BGZF indexes across multiple readers.

## Features

1. **Memory Efficiency**: Both FASTA index and BGZF index are loaded only once and shared.

2. **Thread Safety**: 
   - Read-only shared metadata
   - Per-thread BGZF handles
   - Reference counting with mutex protection

3. **Compatibility**: Follows the same API patterns as the original htslib code.

4. **Independence**: Only needs htslib for specific functions but encapsulates the critical functionality.

5. **Completeness**: Implements all essential functionality from the original faidx API.

## Installation

### Prerequisites
- C compiler (GCC or Clang)
- htslib (installed and available in your system)
- pthread support

### Building and Installing

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/faigz.git
   cd faigz
   ```

2. Build the example and tests:
   ```
   make
   ```

3. Install the header file:
   ```
   sudo make install
   ```
   
   This will install the header to /usr/local/include by default.
   
   To install to a different location:
   ```
   make PREFIX=/your/custom/path install
   ```

4. Run tests with your FASTA file:
   ```
   ./test path/to/your/reference.fa
   ```

## Usage

Include the header in your C/C++ code:

```c
#define REENTRANT_FAIDX_IMPLEMENTATION  // Include this only once in your project
#include <faigz.h>

// Your code here
```

## Usage

### Command Line Tool

The `bench_faigz` tool is provided to measure the performance of the library when accessing BGZF-compressed FASTA files concurrently:

```
Usage: ./bench_faigz [options] <fasta_file>
Options:
  -t INT    Number of threads [4]
  -n INT    Number of sequences to fetch per thread [1000]
  -l INT    Length of each sequence to fetch [100]
  -o FILE   Output fetched sequences to file [none]
  -s INT    Random seed [42]
  -v        Verbose output
  -h        Show this help message
```

Example command:
```bash
# Run with 8 threads, fetching 5000 sequences per thread
./bench_faigz -t 8 -n 5000 -l 200 -v path/to/your/genome.fa.gz
```

### Library Usage

```c
/* Load the index metadata once */
faidx_meta_t *meta = faidx_meta_load(fasta_file, FAI_FASTA, FAI_CREATE);

/* Create a reader for each thread */
faidx_reader_t *reader = faidx_reader_create(meta);

/* Fetch sequence data */
char *seq = faidx_reader_fetch_seq(reader, "chr1", 1000, 1100, &len);

/* Clean up */
faidx_reader_destroy(reader);
faidx_meta_destroy(meta);
```

## API Documentation

### Metadata Functions

- `faidx_meta_t *faidx_meta_load(const char *filename, enum fai_format_options format, int flags)`: Load FASTA/FASTQ index metadata
- `faidx_meta_t *faidx_meta_ref(faidx_meta_t *meta)`: Increment reference count
- `void faidx_meta_destroy(faidx_meta_t *meta)`: Decrement reference count and free if zero
- `int faidx_meta_nseq(const faidx_meta_t *meta)`: Get number of sequences
- `const char *faidx_meta_iseq(const faidx_meta_t *meta, int i)`: Get name of i-th sequence
- `hts_pos_t faidx_meta_seq_len(const faidx_meta_t *meta, const char *seq)`: Get sequence length
- `int faidx_meta_has_seq(const faidx_meta_t *meta, const char *seq)`: Check if sequence exists

### Reader Functions

- `faidx_reader_t *faidx_reader_create(faidx_meta_t *meta)`: Create a reader from shared metadata
- `void faidx_reader_destroy(faidx_reader_t *reader)`: Destroy a reader
- `char *faidx_reader_fetch_seq(faidx_reader_t *reader, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len)`: Fetch sequence
- `char *faidx_reader_fetch_qual(faidx_reader_t *reader, const char *c_name, hts_pos_t p_beg_i, hts_pos_t p_end_i, hts_pos_t *len)`: Fetch quality string (FASTQ only)

## License

This project is licensed under the MIT License - see the LICENSE file for details.
