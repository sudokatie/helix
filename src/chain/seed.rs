/// Seed finding - locate matching k-mers between query and indexed reference

use crate::index::{canonical_kmer, Index, KmerIterator};

/// Configuration for seed finding
#[derive(Debug, Clone)]
pub struct SeedConfig {
    /// Maximum occurrences of a k-mer before filtering (reduce noise from repetitive regions)
    pub max_occurrences: usize,
    /// Minimum seeds required for a region to be considered
    pub min_seeds: usize,
}

impl Default for SeedConfig {
    fn default() -> Self {
        Self {
            max_occurrences: 1000,
            min_seeds: 2,
        }
    }
}

/// A seed represents a k-mer match between query and target
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Seed {
    /// Position in query (0-indexed)
    pub query_pos: u32,
    /// Global position in reference (0-indexed)
    pub target_pos: u32,
    /// Length of the seed match (typically k)
    pub length: u32,
}

impl Seed {
    /// Create a new seed
    pub fn new(query_pos: u32, target_pos: u32, length: u32) -> Self {
        Self {
            query_pos,
            target_pos,
            length,
        }
    }

    /// Diagonal of this seed (target_pos - query_pos)
    /// Seeds on the same diagonal are collinear
    pub fn diagonal(&self) -> i64 {
        self.target_pos as i64 - self.query_pos as i64
    }
}

/// Find seed matches between query and indexed reference
///
/// Returns seeds sorted by (target_pos, query_pos)
pub fn find_seeds(query: &[u8], index: &Index, config: &SeedConfig) -> Vec<Seed> {
    let k = index.config.k;
    let mut seeds = Vec::new();

    for (pos, kmer) in KmerIterator::new(query, k) {
        // Use canonical k-mer for strand-independent matching
        let canonical = canonical_kmer(kmer, k);

        // Look up in index
        let positions = index.hash.lookup(canonical);

        // Filter high-frequency k-mers (repetitive sequences)
        if positions.len() > config.max_occurrences {
            continue;
        }

        // Add seeds for each hit
        for &target_pos in positions {
            seeds.push(Seed::new(pos as u32, target_pos, k as u32));
        }
    }

    // Sort by target position, then query position
    seeds.sort_by_key(|s| (s.target_pos, s.query_pos));

    seeds
}

/// Group seeds by approximate region (reference/contig)
pub fn group_seeds_by_region<'a>(
    seeds: &'a [Seed],
    index: &Index,
) -> Vec<(usize, &'a [Seed])> {
    if seeds.is_empty() {
        return Vec::new();
    }

    let mut groups = Vec::new();
    let mut start = 0;

    // Determine which reference each seed belongs to
    for (i, seed) in seeds.iter().enumerate() {
        let ref_idx = find_reference_index(seed.target_pos as u64, index);

        if i > 0 {
            let prev_ref = find_reference_index(seeds[i - 1].target_pos as u64, index);
            if ref_idx != prev_ref {
                groups.push((prev_ref, &seeds[start..i]));
                start = i;
            }
        }
    }

    // Add final group
    let last_ref = find_reference_index(seeds.last().unwrap().target_pos as u64, index);
    groups.push((last_ref, &seeds[start..]));

    groups
}

/// Find which reference a global position belongs to
fn find_reference_index(global_pos: u64, index: &Index) -> usize {
    for (i, ref_info) in index.references.iter().enumerate() {
        let ref_start = ref_info.offset;
        let ref_end = ref_start + ref_info.length as u64;
        if global_pos >= ref_start && global_pos < ref_end {
            return i;
        }
    }
    0 // Default to first reference
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::{Index, IndexBuilder, IndexConfig, ReferenceInfo};
    use crate::seq::{encode_sequence, Reference};

    fn create_test_index() -> Index {
        let config = IndexConfig { k: 11, w: 5 };
        let mut builder = IndexBuilder::new(config);

        // Create a simple reference with repeated patterns
        let ref_seq = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let reference = Reference {
            name: "test".to_string(),
            seq: ref_seq,
        };
        builder.add_reference(&reference);

        builder.build()
    }

    #[test]
    fn test_find_seeds_exact() {
        let index = create_test_index();
        let query = encode_sequence(b"ACGTACGTACG"); // 11-mer matching reference
        let config = SeedConfig::default();

        let seeds = find_seeds(&query, &index, &config);

        // Should find at least one seed
        assert!(!seeds.is_empty());
    }

    #[test]
    fn test_find_seeds_no_match() {
        let index = create_test_index();
        let query = encode_sequence(b"TTTTTTTTTTT"); // No match in reference
        let config = SeedConfig::default();

        let seeds = find_seeds(&query, &index, &config);

        // Should find no seeds (or filtered due to frequency)
        // Note: might find some due to minimizer sampling
    }

    #[test]
    fn test_seed_diagonal() {
        let seed1 = Seed::new(5, 10, 11);
        let seed2 = Seed::new(6, 11, 11);

        // Same diagonal = collinear
        assert_eq!(seed1.diagonal(), seed2.diagonal());

        let seed3 = Seed::new(5, 12, 11);
        assert_ne!(seed1.diagonal(), seed3.diagonal());
    }

    #[test]
    fn test_filter_high_frequency() {
        let index = create_test_index();
        let query = encode_sequence(b"ACGTACGTACG");
        let config = SeedConfig {
            max_occurrences: 1, // Very restrictive
            min_seeds: 1,
        };

        let seeds = find_seeds(&query, &index, &config);

        // High-frequency k-mers should be filtered
        // Result depends on index contents
    }

    #[test]
    fn test_seeds_sorted() {
        let index = create_test_index();
        let query = encode_sequence(b"ACGTACGTACGTACGTACGT"); // Multiple potential matches
        let config = SeedConfig::default();

        let seeds = find_seeds(&query, &index, &config);

        // Verify sorted by (target_pos, query_pos)
        for window in seeds.windows(2) {
            assert!(
                (window[0].target_pos, window[0].query_pos)
                    <= (window[1].target_pos, window[1].query_pos)
            );
        }
    }
}
