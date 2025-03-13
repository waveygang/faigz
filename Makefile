# faigz - Reentrant FASTA/BGZF index library
# Makefile

# Installation directories
PREFIX = /usr/local
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -O2 -pthread
LDFLAGS = -pthread

# htslib integration
HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null || echo "-I/usr/local/include")
HTSLIB_LIBS := $(shell pkg-config --libs htslib 2>/dev/null || echo "-L/usr/local/lib -lhts")

# Sources and targets
HEADERS = faigz.h
MAIN_SRC = bench_faigz.c
MAIN = faigz

.PHONY: all clean install uninstall

all: $(MAIN)

$(MAIN): $(MAIN_SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(HTSLIB_CFLAGS) -o $@ $< $(HTSLIB_LIBS) $(LDFLAGS)

install: $(HEADERS)
	mkdir -p $(INCLUDEDIR)
	cp $(HEADERS) $(INCLUDEDIR)/
	@echo "Installed faigz.h to $(INCLUDEDIR)"

uninstall:
	rm -f $(INCLUDEDIR)/$(HEADERS)
	@echo "Uninstalled faigz.h from $(INCLUDEDIR)"

clean:
	rm -f $(MAIN) *.o *.gch
