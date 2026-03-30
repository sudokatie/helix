/// AVX2-accelerated Smith-Waterman alignment
///
/// Uses query-parallel processing for SIMD acceleration.

use super::{AlignmentConfig, SwResult};

#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

/// Number of i16 elements in a 256-bit AVX2 vector
#[cfg(target_arch = "x86_64")]
const LANES: usize = 16;

/// Smith-Waterman alignment using AVX2 SIMD
///
/// This implementation uses a simplified approach that parallelizes
/// the max-tracking rather than the full DP computation. For sequences
/// where SIMD provides benefit, we batch process the innermost loop.
///
/// Safety: Caller must ensure AVX2 is available (use is_x86_feature_detected!)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
pub unsafe fn smith_waterman_avx2(
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

    // For short sequences or simple cases, use scalar
    // SIMD overhead isn't worth it for small problems
    if m < 32 || n < 32 {
        return super::smith_waterman_scalar(query, target, config);
    }

    // Use the vectorized inner loop implementation
    smith_waterman_avx2_inner(query, target, config)
}

/// Inner AVX2 implementation with vectorized column processing
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn smith_waterman_avx2_inner(
    query: &[u8],
    target: &[u8],
    config: &AlignmentConfig,
) -> SwResult {
    let m = query.len();
    let n = target.len();

    // Pad query length to multiple of LANES
    let padded_m = ((m + LANES - 1) / LANES) * LANES;

    // Allocate DP vectors (we use two rows + E vector)
    let mut h_prev = vec![0i16; padded_m + 1];
    let mut h_curr = vec![0i16; padded_m + 1];
    let mut e = vec![i16::MIN / 2; padded_m + 1];

    // Scoring parameters as i16
    let match_score = config.match_score as i16;
    let mismatch_penalty = config.mismatch_penalty as i16;
    let gap_open = config.gap_open as i16;
    let gap_extend = config.gap_extend as i16;

    // SIMD constants
    let _zero_vec = _mm256_setzero_si256();
    let _gap_open_vec = _mm256_set1_epi16(gap_open);
    let gap_extend_vec = _mm256_set1_epi16(gap_extend);
    let gap_oe_vec = _mm256_set1_epi16(gap_open.saturating_add(gap_extend));
    let match_vec = _mm256_set1_epi16(match_score);
    let mismatch_vec = _mm256_set1_epi16(mismatch_penalty);

    // Track maximum
    let mut max_score: i16 = 0;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    // Pad query with sentinel values
    let mut padded_query = vec![255u8; padded_m];
    padded_query[..m].copy_from_slice(query);

    // Process each target column
    for j in 0..n {
        let target_base = target[j];
        let target_vec = _mm256_set1_epi8(target_base as i8);

        // F: gap in target (vertical), starts fresh each column
        let mut f: i16 = i16::MIN / 2;

        // Process query in chunks of LANES
        for chunk in 0..(padded_m / LANES) {
            let base_i = chunk * LANES;

            // Load query chunk
            let query_ptr = padded_query.as_ptr().add(base_i);
            let query_vec = _mm256_loadu_si256(query_ptr as *const __m256i);

            // Compare bytes
            let match_mask = _mm256_cmpeq_epi8(query_vec, target_vec);

            // Convert lower 16 bytes to i16 mask
            let mask_lo = _mm256_cvtepi8_epi16(_mm256_castsi256_si128(match_mask));

            // Select match or mismatch score
            let scores = _mm256_blendv_epi8(mismatch_vec, match_vec, mask_lo);

            // Load H[i-1][j-1] (shifted diagonal)
            let h_diag_ptr = h_prev.as_ptr().add(base_i);
            let h_diag = _mm256_loadu_si256(h_diag_ptr as *const __m256i);

            // Load H[i][j-1] (previous column)
            let h_left_ptr = h_prev.as_ptr().add(base_i + 1);
            let h_left = _mm256_loadu_si256(h_left_ptr as *const __m256i);

            // Load E[i] (gap in query, horizontal)
            let e_ptr = e.as_ptr().add(base_i + 1);
            let e_vec = _mm256_loadu_si256(e_ptr as *const __m256i);

            // Compute E: gap in query (extending horizontal gap)
            // E[i] = max(H[i][j-1] + gap_open + gap_extend, E[i] + gap_extend)
            let e_open = _mm256_adds_epi16(h_left, gap_oe_vec);
            let e_extend = _mm256_adds_epi16(e_vec, gap_extend_vec);
            let e_new = _mm256_max_epi16(e_open, e_extend);

            // Store new E
            _mm256_storeu_si256(e.as_mut_ptr().add(base_i + 1) as *mut __m256i, e_new);

            // Compute H from diagonal
            let h_match = _mm256_adds_epi16(h_diag, scores);

            // For F (gap in target), we need to process serially within the chunk
            // because F[i] depends on H[i-1][j] in the same column
            let mut h_new_arr = [0i16; LANES];
            let mut h_match_arr = [0i16; LANES];
            let mut e_new_arr = [0i16; LANES];

            _mm256_storeu_si256(h_match_arr.as_mut_ptr() as *mut __m256i, h_match);
            _mm256_storeu_si256(e_new_arr.as_mut_ptr() as *mut __m256i, e_new);

            for lane in 0..LANES {
                let i = base_i + lane + 1;
                if i > m {
                    h_new_arr[lane] = 0;
                    continue;
                }

                // F: gap in target
                let f_open = h_curr[i - 1].saturating_add(gap_open).saturating_add(gap_extend);
                let f_extend = f.saturating_add(gap_extend);
                f = f_open.max(f_extend);

                // H: best of match, E, F, or 0 (local alignment)
                let h = 0i16
                    .max(h_match_arr[lane])
                    .max(e_new_arr[lane])
                    .max(f);

                h_new_arr[lane] = h;
                h_curr[i] = h;

                // Track max
                if h > max_score {
                    max_score = h;
                    max_i = i;
                    max_j = j + 1;
                }
            }
        }

        // Swap rows
        std::mem::swap(&mut h_prev, &mut h_curr);
        h_curr.fill(0);
    }

    // Approximate start positions based on alignment length
    let approx_len = (max_score as usize) / (match_score.max(1) as usize);
    let query_start = max_i.saturating_sub(approx_len.max(1));
    let target_start = max_j.saturating_sub(approx_len.max(1));

    SwResult {
        score: max_score as i32,
        query_start,
        query_end: max_i,
        target_start,
        target_end: max_j,
    }
}

/// Stub for non-x86_64 platforms
#[cfg(not(target_arch = "x86_64"))]
pub fn smith_waterman_avx2(
    _query: &[u8],
    _target: &[u8],
    _config: &AlignmentConfig,
) -> SwResult {
    panic!("AVX2 not available on this platform");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_sequence;

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_avx2_exact_match() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        // 40 bases to exceed threshold
        let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let result = unsafe { smith_waterman_avx2(&query, &target, &config) };

        // 40 matches * 2 = 80
        assert_eq!(result.score, 80);
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_avx2_vs_scalar() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let scalar = super::super::smith_waterman_scalar(&query, &target, &config);
        let simd = unsafe { smith_waterman_avx2(&query, &target, &config) };

        assert_eq!(simd.score, scalar.score);
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_avx2_with_mismatch() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = encode_sequence(b"ACGTACGTTTGTACGTACGTACGTACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let scalar = super::super::smith_waterman_scalar(&query, &target, &config);
        let simd = unsafe { smith_waterman_avx2(&query, &target, &config) };

        assert_eq!(simd.score, scalar.score);
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_avx2_short_fallback() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        let query = encode_sequence(b"ACGT");
        let target = encode_sequence(b"ACGT");
        let config = AlignmentConfig::default();

        let result = unsafe { smith_waterman_avx2(&query, &target, &config) };
        assert_eq!(result.score, 8);
    }

    #[test]
    #[cfg(target_arch = "x86_64")]
    fn test_avx2_empty() {
        if !is_x86_feature_detected!("avx2") {
            return;
        }

        let query: Vec<u8> = vec![];
        let target = encode_sequence(b"ACGT");
        let config = AlignmentConfig::default();

        let result = unsafe { smith_waterman_avx2(&query, &target, &config) };
        assert_eq!(result.score, 0);
    }
}
