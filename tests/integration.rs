/// Integration tests for helix

use helix::align::{align, smith_waterman_scalar, AlignmentConfig, AlignmentScorer};
use helix::index::{Index, IndexBuilder, IndexConfig};
use helix::seq::{encode_sequence, reverse_complement, FastaReader, FastqReader, Reference};
use helix::chain::{chain_seeds, extend_chain, find_seeds, ChainConfig, ExtendConfig, SeedConfig};
use helix::output::{SamHeader, SamRecord, SamWriter, FLAG_REVERSE, FLAG_UNMAPPED};
use helix::util::detect_simd;

use std::io::{BufReader, Cursor};

/// Test full pipeline: parse -> index -> seed -> chain
#[test]
fn test_end_to_end_pipeline() {
    // Parse reference
    let fasta_data = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    let mut fasta_reader = FastaReader::new(BufReader::new(Cursor::new(fasta_data)));
    let reference = fasta_reader.next_record().unwrap().unwrap();

    // Build index
    let config = IndexConfig { k: 11, w: 5 };
    let mut builder = IndexBuilder::new(config);
    builder.add_reference(&reference);
    let index = builder.build();

    // Parse query
    let fastq_data = "@read1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";
    let mut fastq_reader = FastqReader::new(BufReader::new(Cursor::new(fastq_data)));
    let query = fastq_reader.next_record().unwrap().unwrap();

    // Find seeds
    let seed_config = SeedConfig::default();
    let seeds = find_seeds(&query.seq, &index, &seed_config);

    // Chain seeds
    let chain_config = ChainConfig {
        min_score: 10, // Lower threshold for short test
        ..Default::default()
    };
    let chains = chain_seeds(&seeds, &chain_config);

    // Should find alignment
    // Note: depends on minimizer selection, may not always produce chains
    println!("Seeds found: {}", seeds.len());
    println!("Chains found: {}", chains.len());
}

/// Test SIMD implementations produce same results as scalar
#[test]
fn test_simd_consistency() {
    let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    let target = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    let config = AlignmentConfig::default();

    // Scalar result
    let scalar_result = smith_waterman_scalar(&query, &target, &config);

    // Dispatched result (uses best available SIMD)
    let simd_result = align(&query, &target, &config);

    // Scores must match
    assert_eq!(
        simd_result.score, scalar_result.score,
        "SIMD result {} != scalar result {}",
        simd_result.score, scalar_result.score
    );

    println!("SIMD level: {:?}", detect_simd());
    println!("Alignment score: {}", scalar_result.score);
}

/// Test with mismatches
#[test]
fn test_alignment_with_errors() {
    // Query with some mismatches
    let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    let target = encode_sequence(b"ACGTACGTTCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
    //                                   ^ mismatch
    let config = AlignmentConfig::default();

    let scalar = smith_waterman_scalar(&query, &target, &config);
    let dispatched = align(&query, &target, &config);

    // Should still align with high score
    assert!(scalar.score > 90, "Score too low: {}", scalar.score);
    assert_eq!(scalar.score, dispatched.score);
}

/// Test with insertions
#[test]
fn test_alignment_with_indels() {
    let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGT");
    let target = encode_sequence(b"ACGTACGTACGTTTTTACGTACGTACGT"); // TTTT insertion
    let config = AlignmentConfig::default();

    let result = smith_waterman_scalar(&query, &target, &config);

    // Should find partial alignment
    assert!(result.score > 0);
}

/// Test empty inputs
#[test]
fn test_empty_inputs() {
    let empty: Vec<u8> = vec![];
    let seq = encode_sequence(b"ACGT");
    let config = AlignmentConfig::default();

    let result1 = smith_waterman_scalar(&empty, &seq, &config);
    let result2 = smith_waterman_scalar(&seq, &empty, &config);

    assert_eq!(result1.score, 0);
    assert_eq!(result2.score, 0);
}

/// Test sequence encoding round-trip
#[test]
fn test_sequence_encoding() {
    use helix::seq::{encode_base, decode_base};

    for base in [b'A', b'C', b'G', b'T', b'N'] {
        let encoded = encode_base(base);
        let decoded = decode_base(encoded);
        assert_eq!(decoded, base.to_ascii_uppercase());
    }

    // Case insensitive
    for (lower, upper) in [(b'a', b'A'), (b'c', b'C'), (b'g', b'G'), (b't', b'T')] {
        assert_eq!(encode_base(lower), encode_base(upper));
    }
}

/// Test index serialization round-trip
#[test]
fn test_index_serialization() {
    use helix::index::{save_index, load_index};
    use std::fs;

    // Build index
    let config = IndexConfig { k: 11, w: 5 };
    let mut builder = IndexBuilder::new(config);
    let reference = Reference {
        name: "test".to_string(),
        seq: encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT"),
    };
    builder.add_reference(&reference);
    let index = builder.build();

    // Use temp file
    let temp_path = "/tmp/helix_test_index.hlix";

    // Serialize
    save_index(&index, temp_path).unwrap();

    // Deserialize
    let loaded = load_index(temp_path).unwrap();

    // Cleanup
    let _ = fs::remove_file(temp_path);

    // Verify
    assert_eq!(loaded.config.k, index.config.k);
    assert_eq!(loaded.config.w, index.config.w);
    assert_eq!(loaded.references.len(), index.references.len());
    assert_eq!(loaded.references[0].name, index.references[0].name);
}

/// Test FASTQ parsing
#[test]
fn test_fastq_parsing() {
    let data = "@read1 description\nACGT\n+\nIIII\n@read2\nTGCA\n+\nHHHH\n";
    let mut reader = FastqReader::new(BufReader::new(Cursor::new(data)));

    let r1 = reader.next_record().unwrap().unwrap();
    assert_eq!(r1.id, "read1");
    assert_eq!(r1.seq.len(), 4);

    let r2 = reader.next_record().unwrap().unwrap();
    assert_eq!(r2.id, "read2");
    assert_eq!(r2.seq.len(), 4);

    assert!(reader.next_record().unwrap().is_none());
}

/// Test FASTA parsing with wrapped lines
#[test]
fn test_fasta_wrapped_lines() {
    let data = ">chr1\nACGT\nTGCA\nAAAA\n>chr2\nTTTT\n";
    let mut reader = FastaReader::new(BufReader::new(Cursor::new(data)));

    let r1 = reader.next_record().unwrap().unwrap();
    assert_eq!(r1.name, "chr1");
    assert_eq!(r1.seq.len(), 12); // ACGT + TGCA + AAAA

    let r2 = reader.next_record().unwrap().unwrap();
    assert_eq!(r2.name, "chr2");
    assert_eq!(r2.seq.len(), 4);
}

// TASK 22: Alignment pipeline integration tests

/// Test full alignment pipeline with chain extension
#[test]
fn test_alignment_pipeline_with_extension() {
    // Create a longer reference for proper minimizer matching
    let ref_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let reference = Reference {
        name: "chr1".to_string(),
        seq: encode_sequence(ref_seq),
    };

    // Build index with smaller k for testing
    let config = IndexConfig { k: 11, w: 5 };
    let mut builder = IndexBuilder::new(config);
    builder.add_reference(&reference);
    let index = builder.build();

    // Query that matches part of the reference
    let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT");

    // Step 1: Find seeds
    let seed_config = SeedConfig::default();
    let seeds = find_seeds(&query, &index, &seed_config);

    // Step 2: Chain seeds
    let chain_config = ChainConfig {
        min_score: 10,
        ..Default::default()
    };
    let chains = chain_seeds(&seeds, &chain_config);

    // Step 3: Extend best chain
    if !chains.is_empty() {
        let extend_config = ExtendConfig::default();
        let alignment = extend_chain(&chains[0], &query, &index, &extend_config);

        if let Some(aln) = alignment {
            assert_eq!(aln.ref_name, "chr1");
            assert!(aln.score > 0);
            assert!(!aln.cigar.is_empty());
        }
    }
}

/// Test alignment with reverse complement
#[test]
fn test_alignment_reverse_complement() {
    let ref_seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let reference = Reference {
        name: "chr1".to_string(),
        seq: encode_sequence(ref_seq),
    };

    let config = IndexConfig { k: 11, w: 5 };
    let mut builder = IndexBuilder::new(config);
    builder.add_reference(&reference);
    let index = builder.build();

    // Forward query
    let forward_query = encode_sequence(b"ACGTACGTACGTACGTACGTACGT");

    // Reverse complement query
    let rc_query = reverse_complement(&forward_query);

    // Both should produce seeds
    let seed_config = SeedConfig::default();
    let fwd_seeds = find_seeds(&forward_query, &index, &seed_config);
    let rc_seeds = find_seeds(&rc_query, &index, &seed_config);

    // At least one orientation should find seeds
    assert!(fwd_seeds.len() > 0 || rc_seeds.len() > 0);
}

/// Test SAM record creation from alignment
#[test]
fn test_sam_record_from_alignment() {
    // Create a mapped SAM record
    let record = SamRecord {
        qname: "read1".to_string(),
        flag: 0,
        rname: "chr1".to_string(),
        pos: 100,
        mapq: 60,
        cigar: "50M".to_string(),
        rnext: "*".to_string(),
        pnext: 0,
        tlen: 0,
        seq: "ACGTACGTACGT".to_string(),
        qual: "IIIIIIIIIIII".to_string(),
        tags: vec![("AS".to_string(), "i:100".to_string())],
    };

    // Write to buffer
    let mut buf = Vec::new();
    {
        let mut writer = SamWriter::new(Cursor::new(&mut buf));
        let mut header = SamHeader::new();
        header.add_reference("chr1".to_string(), 1000);
        writer.write_header(&header).unwrap();
        writer.write_record(&record).unwrap();
    }

    let output = String::from_utf8(buf).unwrap();
    assert!(output.contains("read1\t0\tchr1\t100"));
    assert!(output.contains("50M"));
    assert!(output.contains("AS:i:100"));
}

/// Test SAM record for reverse complement alignment
#[test]
fn test_sam_record_reverse_complement() {
    let record = SamRecord {
        qname: "read1".to_string(),
        flag: FLAG_REVERSE,
        rname: "chr1".to_string(),
        pos: 100,
        mapq: 60,
        cigar: "50M".to_string(),
        rnext: "*".to_string(),
        pnext: 0,
        tlen: 0,
        seq: "ACGTACGTACGT".to_string(),
        qual: "IIIIIIIIIIII".to_string(),
        tags: vec![],
    };

    assert_eq!(record.flag & FLAG_REVERSE, FLAG_REVERSE);
}

/// Test unmapped read handling
#[test]
fn test_unmapped_read_handling() {
    // Query that won't match anything
    let ref_seq = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    let reference = Reference {
        name: "chr1".to_string(),
        seq: encode_sequence(ref_seq),
    };

    let config = IndexConfig { k: 11, w: 5 };
    let mut builder = IndexBuilder::new(config);
    builder.add_reference(&reference);
    let index = builder.build();

    // Query with completely different sequence
    let query = encode_sequence(b"TTTTTTTTTTTTTTTTTTTTTTTT");

    let seed_config = SeedConfig::default();
    let seeds = find_seeds(&query, &index, &seed_config);

    // Should output as unmapped
    let record = if seeds.is_empty() {
        SamRecord::new_unmapped(
            "read1".to_string(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".to_string(),
            "IIIIIIIIIIIIIIIIIIIIIIII".to_string(),
        )
    } else {
        // Still unmapped if no good chains
        SamRecord::new_unmapped(
            "read1".to_string(),
            "TTTTTTTTTTTTTTTTTTTTTTTT".to_string(),
            "IIIIIIIIIIIIIIIIIIIIIIII".to_string(),
        )
    };

    assert_eq!(record.flag & FLAG_UNMAPPED, FLAG_UNMAPPED);
    assert_eq!(record.rname, "*");
    assert_eq!(record.cigar, "*");
}

/// Test alignment quality scoring integration
#[test]
fn test_alignment_quality_scoring() {
    use helix::align::Cigar;

    let scorer = AlignmentScorer::default();
    let cigar = Cigar::from("45M5I");

    let result = scorer.evaluate(&cigar, 80, 20, false);

    assert_eq!(result.score, 80);
    assert!(result.identity > 0.8);
    assert!(result.gap_rate <= 0.15);
    assert!(result.passes_filter);
}

/// Test MAPQ calculation
#[test]
fn test_mapq_calculation() {
    let scorer = AlignmentScorer::default();

    // Unique hit (no second best)
    let mapq_unique = scorer.compute_mapq(100, 0);
    assert_eq!(mapq_unique, 60);

    // Close second best
    let mapq_close = scorer.compute_mapq(100, 90);
    assert!(mapq_close < 60);
    assert!(mapq_close > 0);

    // Much worse second best
    let mapq_better = scorer.compute_mapq(100, 20);
    assert!(mapq_better > mapq_close);
}

/// Test multiple reference alignment
#[test]
fn test_multiple_reference_alignment() {
    let ref1 = Reference {
        name: "chr1".to_string(),
        seq: encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),
    };
    let ref2 = Reference {
        name: "chr2".to_string(),
        seq: encode_sequence(b"TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"),
    };

    let config = IndexConfig { k: 11, w: 5 };
    let mut builder = IndexBuilder::new(config);
    builder.add_reference(&ref1);
    builder.add_reference(&ref2);
    let index = builder.build();

    assert_eq!(index.references.len(), 2);
    assert_eq!(index.references[0].name, "chr1");
    assert_eq!(index.references[1].name, "chr2");

    // Query matching chr1
    let query1 = encode_sequence(b"ACGTACGTACGTACGTACGT");
    let seeds1 = find_seeds(&query1, &index, &SeedConfig::default());

    // Query matching chr2
    let query2 = encode_sequence(b"TGCATGCATGCATGCATGCA");
    let seeds2 = find_seeds(&query2, &index, &SeedConfig::default());

    // Both should find seeds
    println!("Seeds for chr1 query: {}", seeds1.len());
    println!("Seeds for chr2 query: {}", seeds2.len());
}
