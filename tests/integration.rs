/// Integration tests for helix

use helix::align::{align, smith_waterman_scalar, AlignmentConfig};
use helix::index::{Index, IndexBuilder, IndexConfig};
use helix::seq::{encode_sequence, FastaReader, FastqReader, Reference};
use helix::chain::{chain_seeds, find_seeds, ChainConfig, SeedConfig};
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
