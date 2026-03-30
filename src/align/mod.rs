mod sw;
mod cigar;
mod sw_avx2;
mod sw_neon;

pub use sw::*;
pub use cigar::*;
pub use sw_avx2::smith_waterman_avx2;
pub use sw_neon::smith_waterman_neon;

use crate::util::{detect_simd, SimdLevel};

/// Perform Smith-Waterman alignment using the best available SIMD implementation
///
/// Automatically dispatches to AVX2, NEON, or scalar based on runtime detection.
pub fn align(query: &[u8], target: &[u8], config: &AlignmentConfig) -> SwResult {
    match detect_simd() {
        SimdLevel::Avx2 => {
            #[cfg(target_arch = "x86_64")]
            {
                // Safety: We just verified AVX2 is available
                unsafe { smith_waterman_avx2(query, target, config) }
            }
            #[cfg(not(target_arch = "x86_64"))]
            {
                smith_waterman_scalar(query, target, config)
            }
        }
        SimdLevel::Neon => {
            #[cfg(target_arch = "aarch64")]
            {
                smith_waterman_neon(query, target, config)
            }
            #[cfg(not(target_arch = "aarch64"))]
            {
                smith_waterman_scalar(query, target, config)
            }
        }
        SimdLevel::None => smith_waterman_scalar(query, target, config),
    }
}

#[cfg(test)]
mod dispatch_tests {
    use super::*;
    use crate::seq::encode_sequence;

    #[test]
    fn test_align_dispatch() {
        // Use long enough sequence to trigger SIMD if available
        let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let result = align(&query, &target, &config);

        // 40 matches * 2 = 80
        assert_eq!(result.score, 80);
    }

    #[test]
    fn test_align_matches_scalar() {
        // Long sequence with some differences
        let query = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let target = encode_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
        let config = AlignmentConfig::default();

        let dispatched = align(&query, &target, &config);
        let scalar = smith_waterman_scalar(&query, &target, &config);

        assert_eq!(dispatched.score, scalar.score);
    }

    #[test]
    fn test_align_short() {
        // Short sequence should use scalar
        let query = encode_sequence(b"ACGTACGT");
        let target = encode_sequence(b"ACGTACGT");
        let config = AlignmentConfig::default();

        let result = align(&query, &target, &config);
        assert_eq!(result.score, 16);
    }
}
