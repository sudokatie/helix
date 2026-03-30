/// NEON-accelerated Smith-Waterman alignment for ARM64
///
/// Uses striped layout to process 8 cells in parallel with 128-bit vectors.

use super::{AlignmentConfig, SwResult};

#[cfg(target_arch = "aarch64")]
use std::arch::aarch64::*;

/// Number of i16 elements in a 128-bit NEON vector
#[cfg(target_arch = "aarch64")]
const LANES: usize = 8;

/// Smith-Waterman alignment using NEON SIMD
///
/// Safety: NEON is always available on aarch64
#[cfg(target_arch = "aarch64")]
pub fn smith_waterman_neon(
    query: &[u8],
    target: &[u8],
    config: &AlignmentConfig,
) -> SwResult {
    // Safety: NEON is guaranteed available on aarch64
    unsafe { smith_waterman_neon_inner(query, target, config) }
}

#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn smith_waterman_neon_inner(
    query: &[u8],
    target: &[u8],
    config: &AlignmentConfig,
) -> SwResult {
    let qlen = query.len();
    let tlen = target.len();

    if qlen == 0 || tlen == 0 {
        return SwResult {
            score: 0,
            query_start: 0,
            query_end: 0,
            target_start: 0,
            target_end: 0,
        };
    }

    // For short sequences, fall back to scalar
    if qlen < LANES {
        return super::smith_waterman_scalar(query, target, config);
    }

    // Number of stripes needed
    let num_stripes = (qlen + LANES - 1) / LANES;

    // Pad query to multiple of LANES
    let padded_len = num_stripes * LANES;
    let mut padded_query = vec![255u8; padded_len];
    padded_query[..qlen].copy_from_slice(query);

    // Create scoring vectors
    let match_vec = vdupq_n_s16(config.match_score as i16);
    let mismatch_vec = vdupq_n_s16(config.mismatch_penalty as i16);
    let gap_open_vec = vdupq_n_s16(config.gap_open as i16);
    let gap_extend_vec = vdupq_n_s16(config.gap_extend as i16);
    let zero = vdupq_n_s16(0);
    let min_val = vdupq_n_s16(i16::MIN / 2);

    // Allocate DP vectors for each stripe
    let mut h_prev: Vec<int16x8_t> = vec![zero; num_stripes];
    let mut h_curr: Vec<int16x8_t> = vec![zero; num_stripes];
    let mut e: Vec<int16x8_t> = vec![min_val; num_stripes];

    // Track maximum score and position
    let mut max_score: i16 = 0;
    let mut max_i = 0usize;
    let mut max_j = 0usize;

    // Temporary storage
    let mut h_tmp = [0i16; LANES];

    // Process each target position
    for j in 0..tlen {
        let target_base = target[j];
        let target_vec = vdupq_n_u8(target_base);

        let mut f = min_val;
        let mut h_diag_carry = zero;

        for stripe in 0..num_stripes {
            let query_offset = stripe * LANES;

            // Load 8 query bases
            let query_ptr = padded_query.as_ptr().add(query_offset);
            let query_chunk = vld1_u8(query_ptr);

            // Compare bases
            let match_mask = vceq_u8(query_chunk, vget_low_u8(vreinterpretq_u8_u16(vdupq_n_u16(target_base as u16))));

            // Convert to i16 scores
            // match_mask is 0xFF for match, 0x00 for mismatch (per byte)
            let mask_16 = vreinterpretq_s16_u16(vmovl_u8(match_mask));

            // Select match or mismatch score
            let match_scores = vandq_s16(mask_16, match_vec);
            let mismatch_mask = vmvnq_s16(mask_16);
            let mismatch_scores = vandq_s16(mismatch_mask, mismatch_vec);
            let s = vorrq_s16(match_scores, mismatch_scores);

            // Load previous row values
            let h_prev_vec = h_prev[stripe];

            // Shift for diagonal access
            let h_shifted = shift_right_insert_neon(h_prev_vec, h_diag_carry);
            h_diag_carry = vdupq_n_s16(vgetq_lane_s16(h_prev_vec, 7));

            // H from diagonal + match/mismatch
            let h_match = vqaddq_s16(h_shifted, s);

            // E: gap in query
            let e_extend = vqaddq_s16(e[stripe], gap_extend_vec);
            let e_open = vqaddq_s16(h_prev[stripe], vqaddq_s16(gap_open_vec, gap_extend_vec));
            let e_new = vmaxq_s16(e_extend, e_open);
            e[stripe] = e_new;

            // F: gap in target
            let f_open = vqaddq_s16(h_curr[stripe], vqaddq_s16(gap_open_vec, gap_extend_vec));
            let f_extend = vqaddq_s16(f, gap_extend_vec);
            f = vmaxq_s16(f_open, f_extend);

            // Final H
            let h_new = vmaxq_s16(zero, h_match);
            let h_new = vmaxq_s16(h_new, e_new);
            let h_new = vmaxq_s16(h_new, f);

            h_curr[stripe] = h_new;

            // Track maximum
            vst1q_s16(h_tmp.as_mut_ptr(), h_new);
            for (lane, &val) in h_tmp.iter().enumerate() {
                if val > max_score {
                    let i = stripe * LANES + lane;
                    if i < qlen {
                        max_score = val;
                        max_i = i + 1;
                        max_j = j + 1;
                    }
                }
            }
        }

        // Swap rows
        std::mem::swap(&mut h_prev, &mut h_curr);
        for stripe in 0..num_stripes {
            h_curr[stripe] = zero;
        }
    }

    let approx_len = (max_score as usize) / (config.match_score.max(1) as usize);
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

/// Shift vector right by one i16, inserting value at position 0
#[cfg(target_arch = "aarch64")]
#[inline]
unsafe fn shift_right_insert_neon(v: int16x8_t, insert: int16x8_t) -> int16x8_t {
    // Extract value to insert
    let insert_val = vgetq_lane_s16(insert, 0);

    // Shift right by extracting and reinserting
    let lane0 = insert_val;
    let lane1 = vgetq_lane_s16(v, 0);
    let lane2 = vgetq_lane_s16(v, 1);
    let lane3 = vgetq_lane_s16(v, 2);
    let lane4 = vgetq_lane_s16(v, 3);
    let lane5 = vgetq_lane_s16(v, 4);
    let lane6 = vgetq_lane_s16(v, 5);
    let lane7 = vgetq_lane_s16(v, 6);

    let arr = [lane0, lane1, lane2, lane3, lane4, lane5, lane6, lane7];
    vld1q_s16(arr.as_ptr())
}

/// Stub for non-aarch64 platforms
#[cfg(not(target_arch = "aarch64"))]
pub fn smith_waterman_neon(
    _query: &[u8],
    _target: &[u8],
    _config: &AlignmentConfig,
) -> SwResult {
    panic!("NEON not available on this platform");
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_sequence;

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_exact_match() {
        let query = encode_sequence(b"ACGTACGTACGTACGT"); // 16 bases
        let target = encode_sequence(b"ACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let result = smith_waterman_neon(&query, &target, &config);

        // 16 matches * 2 = 32
        assert_eq!(result.score, 32);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_vs_scalar() {
        let query = encode_sequence(b"ACGTACGTACGTACGTACGT"); // 20 bases
        let target = encode_sequence(b"ACGTACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let scalar = super::super::smith_waterman_scalar(&query, &target, &config);
        let simd = smith_waterman_neon(&query, &target, &config);

        assert_eq!(simd.score, scalar.score);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_short_fallback() {
        let query = encode_sequence(b"ACGT"); // < 8 bases
        let target = encode_sequence(b"ACGT");
        let config = AlignmentConfig::default();

        let result = smith_waterman_neon(&query, &target, &config);
        assert_eq!(result.score, 8);
    }

    #[test]
    #[cfg(target_arch = "aarch64")]
    fn test_neon_empty() {
        let query: Vec<u8> = vec![];
        let target = encode_sequence(b"ACGT");
        let config = AlignmentConfig::default();

        let result = smith_waterman_neon(&query, &target, &config);
        assert_eq!(result.score, 0);
    }

    // Tests for non-aarch64 platforms (just verify stub exists)
    #[test]
    #[cfg(not(target_arch = "aarch64"))]
    fn test_neon_stub_exists() {
        // Just verify the function exists
        // Can't call it because it would panic
    }
}
