# helix

DNA sequence aligner with SIMD acceleration. Because minimap2 is great, but sometimes you want to understand what's actually happening.

## Why This Exists?

Every bioinformatics lab has that one script that calls BWA or Bowtie and hopes for the best. Helix is for when you want to:

- Learn how sequence alignment actually works (Smith-Waterman, seed-and-extend, the whole deal)
- Have a codebase small enough to read in an afternoon
- Get SIMD acceleration without needing a PhD to understand the implementation

It's not trying to replace production aligners. It's trying to be the aligner you can actually modify.

## Features

- FASTQ/FASTA parsing with gzip support
- K-mer indexing with minimizer sampling (because indexing every k-mer is for masochists)
- Smith-Waterman alignment with affine gap penalties
- AVX2/NEON SIMD acceleration (automatic detection)
- Seed-and-extend alignment pipeline
- SAM output compatible with downstream tools

## Quick Start

```bash
# Build
cargo build --release

# Index a reference genome
helix index reference.fa -o ref_index

# Align reads
helix align ref_index reads.fq -o alignments.sam
```

## Installation

```bash
# From source
git clone https://github.com/sudokatie/helix
cd helix
cargo build --release

# Binary will be at target/release/helix
```

Requires Rust 1.70+ (we use edition 2021 features).

## Usage

### Indexing

```bash
helix index <reference.fa> -o <output_prefix> [-k <kmer_size>]
```

Options:
- `-k, --kmer`: K-mer size for indexing (default: 15)
- `-o, --output`: Output prefix for index files

Creates `<prefix>.hlix` index file.

### Alignment

```bash
helix align <index_prefix> <reads.fq> -o <output.sam> [-t <threads>]
```

Options:
- `-o, --output`: Output SAM file
- `-t, --threads`: Number of threads (default: 1)

### Info

```bash
helix info <index_prefix>
```

Shows index statistics (reference count, k-mer size, etc).

## Algorithm

Helix uses a seed-and-extend approach:

1. **Indexing**: Extract minimizers from reference, build hash table
2. **Seeding**: Find matching k-mers between query and index
3. **Chaining**: Connect collinear seeds using dynamic programming
4. **Extension**: Extend chains with banded Smith-Waterman

The Smith-Waterman implementation supports affine gap penalties and has SIMD-accelerated versions for:
- x86-64 with AVX2 (256-bit, 16 cells parallel)
- ARM64 with NEON (128-bit, 8 cells parallel)

Falls back to scalar automatically when SIMD isn't available.

## Performance

On a modern x86-64 machine with AVX2:
- Indexing: ~10 MB/s reference sequence
- Alignment: depends heavily on read length and error rate

For production workloads, use minimap2. For learning and experimentation, helix is more approachable.

## Project Structure

```
src/
  seq/       - Sequence encoding, FASTQ/FASTA parsing
  index/     - K-mer extraction, minimizers, hash index
  align/     - Smith-Waterman (scalar + SIMD), CIGAR
  chain/     - Seed finding, chaining, extension
  output/    - SAM format output
  util/      - SIMD detection
```

## License

MIT

## Author

Katie

---

*Built to understand, not to compete.*
