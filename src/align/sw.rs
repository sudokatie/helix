/// Configuration for alignment scoring
#[derive(Debug, Clone)]
pub struct AlignmentConfig {
    /// Score for matching bases (positive)
    pub match_score: i32,
    /// Penalty for mismatching bases (negative or zero)
    pub mismatch_penalty: i32,
    /// Penalty for opening a gap (negative)
    pub gap_open: i32,
    /// Penalty for extending a gap (negative)
    pub gap_extend: i32,
}

impl Default for AlignmentConfig {
    fn default() -> Self {
        Self {
            match_score: 2,
            mismatch_penalty: -1,
            gap_open: -2,
            gap_extend: -1,
        }
    }
}

/// Result of Smith-Waterman alignment
#[derive(Debug, Clone)]
pub struct SwResult {
    /// Alignment score
    pub score: i32,
    /// Start position in query (0-indexed)
    pub query_start: usize,
    /// End position in query (exclusive)
    pub query_end: usize,
    /// Start position in target (0-indexed)
    pub target_start: usize,
    /// End position in target (exclusive)
    pub target_end: usize,
}

/// Scalar Smith-Waterman local alignment with affine gaps
pub fn smith_waterman_scalar(
    query: &[u8],
    target: &[u8],
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

    // DP matrices
    // H[i][j] = best score ending at (i,j) with match/mismatch
    // E[i][j] = best score ending at (i,j) with gap in query (deletion)
    // F[i][j] = best score ending at (i,j) with gap in target (insertion)
    let mut h = vec![vec![0i32; n + 1]; m + 1];
    let mut e = vec![vec![i32::MIN / 2; n + 1]; m + 1];
    let mut f = vec![vec![i32::MIN / 2; n + 1]; m + 1];

    let mut max_score = 0;
    let mut max_i = 0;
    let mut max_j = 0;

    for i in 1..=m {
        for j in 1..=n {
            // Match/mismatch score
            let s = if query[i - 1] == target[j - 1] {
                config.match_score
            } else {
                config.mismatch_penalty
            };

            // E: gap in query (horizontal move in matrix)
            e[i][j] = (h[i][j - 1] + config.gap_open + config.gap_extend)
                .max(e[i][j - 1] + config.gap_extend);

            // F: gap in target (vertical move in matrix)
            f[i][j] = (h[i - 1][j] + config.gap_open + config.gap_extend)
                .max(f[i - 1][j] + config.gap_extend);

            // H: best score at this cell
            h[i][j] = 0
                .max(h[i - 1][j - 1] + s)
                .max(e[i][j])
                .max(f[i][j]);

            // Track maximum
            if h[i][j] > max_score {
                max_score = h[i][j];
                max_i = i;
                max_j = j;
            }
        }
    }

    // Traceback to find alignment start
    let (query_start, target_start) = if max_score > 0 {
        traceback(&h, &e, &f, query, target, config, max_i, max_j)
    } else {
        (0, 0)
    };

    SwResult {
        score: max_score,
        query_start,
        query_end: max_i,
        target_start,
        target_end: max_j,
    }
}

/// Traceback to find alignment start position
fn traceback(
    h: &[Vec<i32>],
    e: &[Vec<i32>],
    f: &[Vec<i32>],
    query: &[u8],
    target: &[u8],
    config: &AlignmentConfig,
    end_i: usize,
    end_j: usize,
) -> (usize, usize) {
    let mut i = end_i;
    let mut j = end_j;

    #[derive(Clone, Copy, PartialEq)]
    enum State {
        Match,
        InsertQuery, // gap in query = E
        InsertTarget, // gap in target = F
    }

    let mut state = State::Match;

    while i > 0 && j > 0 && h[i][j] > 0 {
        match state {
            State::Match => {
                let s = if query[i - 1] == target[j - 1] {
                    config.match_score
                } else {
                    config.mismatch_penalty
                };

                if h[i][j] == h[i - 1][j - 1] + s {
                    i -= 1;
                    j -= 1;
                } else if h[i][j] == e[i][j] {
                    state = State::InsertQuery;
                } else if h[i][j] == f[i][j] {
                    state = State::InsertTarget;
                } else {
                    break; // Local alignment ended
                }
            }
            State::InsertQuery => {
                if e[i][j] == h[i][j - 1] + config.gap_open + config.gap_extend {
                    state = State::Match;
                }
                j -= 1;
            }
            State::InsertTarget => {
                if f[i][j] == h[i - 1][j] + config.gap_open + config.gap_extend {
                    state = State::Match;
                }
                i -= 1;
            }
        }
    }

    (i, j)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_sequence;

    #[test]
    fn test_exact_match() {
        let query = encode_sequence(b"ACGT");
        let target = encode_sequence(b"ACGT");
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);

        assert_eq!(result.score, 8); // 4 matches * 2
        assert_eq!(result.query_start, 0);
        assert_eq!(result.query_end, 4);
    }

    #[test]
    fn test_single_mismatch() {
        let query = encode_sequence(b"ACGT");
        let target = encode_sequence(b"ACTT"); // G->T mismatch
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);

        // 3 matches (6) + 1 mismatch (-1) = 5
        assert_eq!(result.score, 5);
    }

    #[test]
    fn test_single_gap() {
        let query = encode_sequence(b"ACGT");
        let target = encode_sequence(b"ACGGT"); // extra G
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);

        // Should find alignment with gap
        assert!(result.score > 0);
    }

    #[test]
    fn test_no_alignment() {
        let query = encode_sequence(b"AAAA");
        let target = encode_sequence(b"TTTT");
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);

        // All mismatches, local alignment score is 0
        assert_eq!(result.score, 0);
    }

    #[test]
    fn test_partial_match() {
        let query = encode_sequence(b"TTACGTAA");
        let target = encode_sequence(b"NNACGTNN");
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);

        // Should find the ACGT match in the middle
        assert!(result.score >= 8); // At least 4 matches
    }

    #[test]
    fn test_empty_sequences() {
        let query: Vec<u8> = vec![];
        let target = encode_sequence(b"ACGT");
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);
        assert_eq!(result.score, 0);

        let result = smith_waterman_scalar(&target, &query, &config);
        assert_eq!(result.score, 0);
    }

    #[test]
    fn test_affine_gap() {
        let query = encode_sequence(b"ACGTACGT");
        let target = encode_sequence(b"ACGTACGT"); // same as query
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);
        assert_eq!(result.score, 16); // 8 matches * 2
    }

    #[test]
    fn test_gap_extension_cheaper() {
        let query = encode_sequence(b"ACGTTTTTACGT");
        let target = encode_sequence(b"ACGTACGT"); // missing TTTTT
        let config = AlignmentConfig {
            match_score: 2,
            mismatch_penalty: -1,
            gap_open: -5,
            gap_extend: -1,
        };

        let result = smith_waterman_scalar(&query, &target, &config);

        // Should prefer one long gap over multiple short gaps
        assert!(result.score > 0);
    }

    #[test]
    fn test_alignment_positions() {
        // Use sequences where the alignment is clearly in the middle
        let query = encode_sequence(b"TTTTACGTAAAA");
        let target = encode_sequence(b"GGGGACGTCCCC");
        let config = AlignmentConfig::default();

        let result = smith_waterman_scalar(&query, &target, &config);

        // The best alignment should be ACGT in the middle
        assert_eq!(result.score, 8); // 4 matches * 2
        // Alignment should start after the mismatching prefix
        assert!(result.query_start >= 4);
        assert!(result.target_start >= 4);
    }
}
