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
EXAMPLE_SRC = example.c
EXAMPLE = example
TEST_SRC = test.c
TEST = test

.PHONY: all clean install test uninstall

all: $(EXAMPLE)

$(EXAMPLE): $(EXAMPLE_SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(HTSLIB_CFLAGS) -o $@ $< $(HTSLIB_LIBS) $(LDFLAGS)

test: $(TEST)
	./$(TEST)

$(TEST): $(TEST_SRC) $(HEADERS)
	$(CC) $(CFLAGS) $(HTSLIB_CFLAGS) -o $@ $< $(HTSLIB_LIBS) $(LDFLAGS)

install: $(HEADERS)
	mkdir -p $(INCLUDEDIR)
	cp $(HEADERS) $(INCLUDEDIR)/
	@echo "Installed faigz.h to $(INCLUDEDIR)"

uninstall:
	rm -f $(INCLUDEDIR)/$(HEADERS)
	@echo "Uninstalled faigz.h from $(INCLUDEDIR)"

clean:
	rm -f $(EXAMPLE) $(TEST) *.o *.gch
