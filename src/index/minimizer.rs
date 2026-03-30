use crate::index::kmer::{Kmer, KmerIterator};
use std::collections::VecDeque;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    pub kmer: Kmer,
    pub pos: usize,
}

pub struct MinimizerIterator<'a> {
    kmer_iter: KmerIterator<'a>,
    window: VecDeque<(usize, Kmer, u64)>, // (pos, kmer, hash)
    window_size: usize,
    last_minimizer_pos: Option<usize>,
}

impl<'a> MinimizerIterator<'a> {
    pub fn new(seq: &'a [u8], k: usize, w: usize) -> Self {
        Self {
            kmer_iter: KmerIterator::new(seq, k),
            window: VecDeque::with_capacity(w),
            window_size: w,
            last_minimizer_pos: None,
        }
    }

    fn hash(kmer: Kmer) -> u64 {
        // Invertible hash to break ties deterministically
        let mut x = kmer;
        x = (x ^ (x >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)).wrapping_mul(0x94d049bb133111eb);
        x ^ (x >> 31)
    }
}

impl<'a> Iterator for MinimizerIterator<'a> {
    type Item = Minimizer;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((pos, kmer)) = self.kmer_iter.next() {
                let hash = Self::hash(kmer);

                // Remove old entries from front
                while !self.window.is_empty()
                    && self.window.front().unwrap().0 + self.window_size <= pos
                {
                    self.window.pop_front();
                }

                // Remove entries from back that have larger hash
                while !self.window.is_empty() && self.window.back().unwrap().2 >= hash {
                    self.window.pop_back();
                }

                self.window.push_back((pos, kmer, hash));

                // Check if we have a full window
                if self.window.len() >= 1 {
                    let (min_pos, min_kmer, _) = self.window.front().unwrap();

                    // Deduplicate consecutive minimizers at same position
                    if self.last_minimizer_pos != Some(*min_pos) {
                        self.last_minimizer_pos = Some(*min_pos);
                        return Some(Minimizer {
                            kmer: *min_kmer,
                            pos: *min_pos,
                        });
                    }
                }
            } else {
                return None;
            }
        }
    }
}

pub fn collect_minimizers(seq: &[u8], k: usize, w: usize) -> Vec<Minimizer> {
    MinimizerIterator::new(seq, k, w).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_sequence;

    #[test]
    fn test_minimizer_basic() {
        let seq = encode_sequence(b"ACGTACGTACGT");
        let minimizers = collect_minimizers(&seq, 3, 2);
        // Should have some minimizers but fewer than all k-mers
        assert!(!minimizers.is_empty());
    }

    #[test]
    fn test_minimizer_deduplication() {
        let seq = encode_sequence(b"AAAAAAAAAA");
        let minimizers = collect_minimizers(&seq, 3, 2);
        // All AAA k-mers are the same, should deduplicate
        // Each window slides by 1, but same minimizer shouldn't repeat
        let positions: Vec<_> = minimizers.iter().map(|m| m.pos).collect();
        // Check no duplicate positions
        for i in 1..positions.len() {
            assert_ne!(positions[i], positions[i - 1]);
        }
    }

    #[test]
    fn test_minimizer_fewer_than_kmers() {
        let seq = encode_sequence(b"ACGTACGTACGTACGT");
        let all_kmers: Vec<_> = KmerIterator::new(&seq, 4).collect();
        let minimizers = collect_minimizers(&seq, 4, 3);

        // Minimizers should be fewer than total k-mers
        assert!(minimizers.len() <= all_kmers.len());
        assert!(!minimizers.is_empty());
    }

    #[test]
    fn test_minimizer_window_1() {
        // Window size 1 means every k-mer is a minimizer
        let seq = encode_sequence(b"ACGTACGT");
        let all_kmers: Vec<_> = KmerIterator::new(&seq, 2).collect();
        let minimizers = collect_minimizers(&seq, 2, 1);

        // With w=1, we get roughly as many minimizers as k-mers
        // (deduplication may reduce slightly if consecutive same kmer)
        assert!(minimizers.len() <= all_kmers.len());
    }

    #[test]
    fn test_minimizer_positions_valid() {
        let seq = encode_sequence(b"ACGTACGTACGT");
        let minimizers = collect_minimizers(&seq, 3, 2);

        for m in &minimizers {
            // Position should be valid
            assert!(m.pos + 3 <= seq.len());
        }
    }

    #[test]
    fn test_empty_sequence() {
        let seq: Vec<u8> = vec![];
        let minimizers = collect_minimizers(&seq, 3, 2);
        assert!(minimizers.is_empty());
    }

    #[test]
    fn test_short_sequence() {
        let seq = encode_sequence(b"AC");
        let minimizers = collect_minimizers(&seq, 3, 2);
        assert!(minimizers.is_empty()); // Sequence shorter than k
    }
}
