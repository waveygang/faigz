// test_faigz.cpp - C++ test for faigz library
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

// Include the faigz.h header with implementation
#define REENTRANT_FAIDX_IMPLEMENTATION
#include "faigz.h"

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <fasta_file>\n";
        return 1;
    }

    const char *fasta_file = argv[1];
    std::cout << "Testing faigz C++ integration with file: " << fasta_file << std::endl;
    
    // Load metadata
    faidx_meta_t *meta = faidx_meta_load(fasta_file, FAI_FASTA, FAI_CREATE);
    if (!meta) {
        std::cerr << "Failed to load FASTA index" << std::endl;
        return 1;
    }
    
    std::cout << "Loaded index with " << faidx_meta_nseq(meta) << " sequences" << std::endl;
    
    // Create a reader
    faidx_reader_t *reader = faidx_reader_create(meta);
    if (!reader) {
        std::cerr << "Failed to create reader" << std::endl;
        faidx_meta_destroy(meta);
        return 1;
    }
    
    // Print some sequence names as a test
    std::cout << "First 5 sequence names (or fewer if less available):" << std::endl;
    int max_to_show = std::min(5, faidx_meta_nseq(meta));
    for (int i = 0; i < max_to_show; i++) {
        const char *seq_name = faidx_meta_iseq(meta, i);
        hts_pos_t seq_len = faidx_meta_seq_len(meta, seq_name);
        std::cout << "  " << i+1 << ": " << seq_name << " (length: " << seq_len << ")" << std::endl;
    }
    
    // Fetch a sequence as a test if sequences are available
    if (faidx_meta_nseq(meta) > 0) {
        const char *seq_name = faidx_meta_iseq(meta, 0);
        std::cout << "\nFetching the first 10 bases from " << seq_name << ":" << std::endl;
        
        hts_pos_t len;
        char *seq = faidx_reader_fetch_seq(reader, seq_name, 0, 9, &len);
        
        if (seq) {
            std::cout << "Sequence: " << seq << " (length: " << len << ")" << std::endl;
            free(seq);
        } else {
            std::cout << "Failed to fetch sequence" << std::endl;
        }
    }
    
    // Test reference counting
    std::cout << "\nTesting reference counting:" << std::endl;
    std::cout << "Creating 5 additional readers..." << std::endl;
    
    // Create a vector to hold readers
    std::vector<faidx_reader_t*> readers;
    
    for (int i = 0; i < 5; i++) {
        faidx_reader_t *r = faidx_reader_create(meta);
        if (r) {
            readers.push_back(r);
            std::cout << "  Created reader " << i+1 << std::endl;
        }
    }
    
    // Destroy the readers
    std::cout << "Destroying readers..." << std::endl;
    for (auto r : readers) {
        faidx_reader_destroy(r);
    }
    readers.clear();
    
    // Clean up
    faidx_reader_destroy(reader);
    faidx_meta_destroy(meta);
    
    std::cout << "\nC++ test completed successfully!" << std::endl;
    return 0;
}
