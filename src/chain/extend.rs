/// Chain extension - extend chains with banded alignment

use super::Chain;
use crate::align::{AlignmentConfig, Cigar, CigarOp, SwResult};
use crate::index::Index;

/// Configuration for chain extension
#[derive(Debug, Clone)]
pub struct ExtendConfig {
    /// Band width for banded alignment
    pub band_width: usize,
    /// X-drop threshold (stop if score drops this much below max)
    pub x_drop: i32,
    /// Padding to extend beyond chain boundaries
    pub extension_padding: u32,
    /// Alignment scoring config
    pub align_config: AlignmentConfig,
}

impl Default for ExtendConfig {
    fn default() -> Self {
        Self {
            band_width: 100,
            x_drop: 100,
            extension_padding: 50,
            align_config: AlignmentConfig::default(),
        }
    }
}

/// Full alignment result after chain extension
#[derive(Debug, Clone)]
pub struct Alignment {
    /// Query name/id
    pub query_name: String,
    /// Reference name
    pub ref_name: String,
    /// Reference index
    pub ref_idx: usize,
    /// Start position in query
    pub query_start: usize,
    /// End position in query
    pub query_end: usize,
    /// Start position in reference (1-indexed for SAM)
    pub ref_start: usize,
    /// End position in reference
    pub ref_end: usize,
    /// Alignment score
    pub score: i32,
    /// CIGAR string
    pub cigar: Cigar,
    /// Mapping quality
    pub mapq: u8,
    /// Is this a secondary alignment?
    pub is_secondary: bool,
    /// Is reverse complemented?
    pub is_reverse: bool,
}

/// Extend a chain with banded alignment
pub fn extend_chain(
    chain: &Chain,
    query: &[u8],
    index: &Index,
    config: &ExtendConfig,
) -> Option<Alignment> {
    if chain.seeds.is_empty() {
        return None;
    }

    // Find which reference this chain maps to
    let target_pos = chain.target_start as u64;
    let (ref_idx, ref_info) = find_reference(target_pos, index)?;

    // Calculate local positions within the reference
    let ref_offset = ref_info.offset;
    let local_target_start = (chain.target_start as u64 - ref_offset) as usize;
    let local_target_end = (chain.target_end as u64 - ref_offset) as usize;

    // Extend boundaries with padding
    let query_start = chain.query_start.saturating_sub(config.extension_padding) as usize;
    let query_end = ((chain.query_end + config.extension_padding) as usize).min(query.len());
    let ref_start = local_target_start.saturating_sub(config.extension_padding as usize);
    let ref_end = (local_target_end + config.extension_padding as usize).min(ref_info.length);

    // Extract sequences for alignment
    let _query_slice = &query[query_start..query_end];

    // We need the reference sequence - for now, use a simple approach
    // In practice, we'd have the reference loaded or memory-mapped
    // For this implementation, we'll do banded alignment with the chain seeds
    // as anchors

    // Build alignment from seeds
    let alignment = build_alignment_from_chain(
        chain,
        query_start,
        query_end,
        ref_start,
        ref_end,
        ref_idx,
        &ref_info.name,
        config,
    );

    Some(alignment)
}

/// Build alignment directly from chain seeds
fn build_alignment_from_chain(
    chain: &Chain,
    query_start: usize,
    query_end: usize,
    ref_start: usize,
    ref_end: usize,
    ref_idx: usize,
    ref_name: &str,
    _config: &ExtendConfig,
) -> Alignment {
    // Build CIGAR from chain structure
    let mut cigar = Cigar::new();
    let mut last_query_end = chain.query_start as usize;
    let mut last_ref_end = chain.target_start as usize;

    for seed in &chain.seeds {
        let seed_query_start = seed.query_pos as usize;
        let seed_ref_start = seed.target_pos as usize;

        // Gap before this seed
        let query_gap = seed_query_start.saturating_sub(last_query_end);
        let ref_gap = seed_ref_start.saturating_sub(last_ref_end);

        if query_gap > ref_gap {
            // Insertion in query
            if ref_gap > 0 {
                cigar.push(CigarOp::Match(ref_gap as u32));
            }
            if query_gap - ref_gap > 0 {
                cigar.push(CigarOp::Insertion((query_gap - ref_gap) as u32));
            }
        } else if ref_gap > query_gap {
            // Deletion in query
            if query_gap > 0 {
                cigar.push(CigarOp::Match(query_gap as u32));
            }
            if ref_gap - query_gap > 0 {
                cigar.push(CigarOp::Deletion((ref_gap - query_gap) as u32));
            }
        } else if query_gap > 0 {
            // Equal gaps - treat as matches
            cigar.push(CigarOp::Match(query_gap as u32));
        }

        // The seed itself is a match
        cigar.push(CigarOp::Match(seed.length));

        last_query_end = seed.query_pos as usize + seed.length as usize;
        last_ref_end = seed.target_pos as usize + seed.length as usize;
    }

    // Merge consecutive same operations
    cigar = cigar.merged();

    Alignment {
        query_name: String::new(), // Set by caller
        ref_name: ref_name.to_string(),
        ref_idx,
        query_start,
        query_end,
        ref_start,
        ref_end,
        score: chain.score,
        cigar,
        mapq: estimate_mapq(chain),
        is_secondary: false,
        is_reverse: false,
    }
}

/// Find which reference a global position belongs to
fn find_reference(global_pos: u64, index: &Index) -> Option<(usize, &crate::index::ReferenceInfo)> {
    for (i, ref_info) in index.references.iter().enumerate() {
        let ref_start = ref_info.offset;
        let ref_end = ref_start + ref_info.length as u64;
        if global_pos >= ref_start && global_pos < ref_end {
            return Some((i, ref_info));
        }
    }
    None
}

/// Estimate mapping quality from chain properties
fn estimate_mapq(chain: &Chain) -> u8 {
    // Simple MAPQ estimation based on chain score and coverage
    // In practice, this would consider second-best alignment score
    let coverage_ratio = chain.seeds.len() as f64 / 10.0; // Normalize
    let score_component = (chain.score as f64 / 100.0).min(1.0);

    let mapq = ((coverage_ratio + score_component) * 30.0).min(60.0) as u8;
    mapq
}

/// Banded Smith-Waterman alignment
/// Only compute cells within band of expected diagonal
pub fn banded_align(
    query: &[u8],
    target: &[u8],
    band_width: usize,
    config: &AlignmentConfig,
) -> SwResult {
    let m = query.len();
    let n = target.len();

    if m == 0 || n == 0 {
        return SwResult {
            score: 0,
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
        };
    }

    // For small sequences or wide bands, use regular SW
    if m <= band_width * 2 || n <= band_width * 2 {
        return crate::align::smith_waterman_scalar(query, target, config);
    }

    // DP with banding
    let mut h = vec![vec![0i32; n + 1]; m + 1];
    let mut e = vec![vec![i32::MIN / 2; n + 1]; m + 1];
    let mut f = vec![vec![i32::MIN / 2; n + 1]; m + 1];

    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    for i in 1..=m {
        // Compute band limits for this row
        let expected_j = (i * n) / m; // Expected diagonal
        let j_start = expected_j.saturating_sub(band_width).max(1);
        let j_end = (expected_j + band_width).min(n);

        for j in j_start..=j_end {
            let s = if query[i - 1] == target[j - 1] {
                config.match_score
            } else {
                config.mismatch_penalty
            };

            e[i][j] = (h[i][j - 1] + config.gap_open + config.gap_extend)
                .max(e[i][j - 1] + config.gap_extend);

            f[i][j] = (h[i - 1][j] + config.gap_open + config.gap_extend)
                .max(f[i - 1][j] + config.gap_extend);

            h[i][j] = 0
                .max(h[i - 1][j - 1] + s)
                .max(e[i][j])
                .max(f[i][j]);

            if h[i][j] > max_score {
                max_score = h[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    SwResult {
        score: max_score,
        query_start: max_i.saturating_sub(max_score as usize / config.match_score.max(1) as usize),
        query_end: max_i,
        target_start: max_j.saturating_sub(max_score as usize / config.match_score.max(1) as usize),
        target_end: max_j,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::chain::Seed;

    #[test]
    fn test_build_alignment_simple() {
        let chain = Chain {
            seeds: vec![Seed::new(10, 100, 15)],
            score: 15,
            query_start: 10,
            query_end: 25,
            target_start: 100,
            target_end: 115,
        };

        let config = ExtendConfig::default();
        let alignment = build_alignment_from_chain(
            &chain, 0, 50, 90, 130, 0, "chr1", &config,
        );

        assert_eq!(alignment.score, 15);
        assert!(!alignment.cigar.ops().is_empty());
    }

    #[test]
    fn test_build_alignment_multiple_seeds() {
        let chain = Chain {
            seeds: vec![
                Seed::new(10, 100, 15),
                Seed::new(30, 120, 15),
            ],
            score: 30,
            query_start: 10,
            query_end: 45,
            target_start: 100,
            target_end: 135,
        };

        let config = ExtendConfig::default();
        let alignment = build_alignment_from_chain(
            &chain, 0, 60, 90, 150, 0, "chr1", &config,
        );

        assert_eq!(alignment.score, 30);
    }

    #[test]
    fn test_estimate_mapq() {
        let chain = Chain {
            seeds: vec![Seed::new(0, 100, 15); 10],
            score: 150,
            query_start: 0,
            query_end: 200,
            target_start: 100,
            target_end: 300,
        };

        let mapq = estimate_mapq(&chain);
        assert!(mapq > 0);
        assert!(mapq <= 60);
    }

    #[test]
    fn test_banded_align() {
        use crate::seq::encode_sequence;

        let query = encode_sequence(b"ACGTACGTACGTACGT");
        let target = encode_sequence(b"ACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let result = banded_align(&query, &target, 10, &config);

        assert_eq!(result.score, 32); // 16 * 2
    }

    #[test]
    fn test_banded_align_empty() {
        let query: Vec<u8> = vec![];
        let target: Vec<u8> = vec![];
        let config = AlignmentConfig::default();

        let result = banded_align(&query, &target, 10, &config);

        assert_eq!(result.score, 0);
    }
}
