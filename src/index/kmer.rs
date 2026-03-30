use crate::seq::{complement, N};

pub type Kmer = u64;

pub struct KmerIterator<'a> {
    seq: &'a [u8],
    k: usize,
    pos: usize,
    current: Kmer,
    valid_bases: usize,
    mask: Kmer,
}

impl<'a> KmerIterator<'a> {
    pub fn new(seq: &'a [u8], k: usize) -> Self {
        assert!(k > 0 && k <= 32, "k must be between 1 and 32");
        Self {
            seq,
            k,
            pos: 0,
            current: 0,
            valid_bases: 0,
            mask: (1u64 << (2 * k)) - 1,
        }
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = (usize, Kmer);

    fn next(&mut self) -> Option<Self::Item> {
        while self.pos < self.seq.len() {
            let base = self.seq[self.pos];
            self.pos += 1;

            if base == N {
                self.valid_bases = 0;
                self.current = 0;
                continue;
            }

            self.current = ((self.current << 2) | (base as u64)) & self.mask;
            self.valid_bases += 1;

            if self.valid_bases >= self.k {
                return Some((self.pos - self.k, self.current));
            }
        }
        None
    }
}

pub fn reverse_complement_kmer(kmer: Kmer, k: usize) -> Kmer {
    let mut rc = 0u64;
    let mut remaining = kmer;
    for _ in 0..k {
        let base = (remaining & 3) as u8;
        rc = (rc << 2) | (complement(base) as u64);
        remaining >>= 2;
    }
    rc
}

pub fn canonical_kmer(kmer: Kmer, k: usize) -> Kmer {
    let rc = reverse_complement_kmer(kmer, k);
    kmer.min(rc)
}

pub fn kmer_to_string(kmer: Kmer, k: usize) -> String {
    let bases = ['A', 'C', 'G', 'T'];
    let mut s = String::with_capacity(k);
    for i in 0..k {
        s.push(bases[((kmer >> (2 * (k - 1 - i))) & 3) as usize]);
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_sequence;

    #[test]
    fn test_kmer_extraction_basic() {
        let seq = encode_sequence(b"ACGT");
        let kmers: Vec<_> = KmerIterator::new(&seq, 2).collect();
        // AC, CG, GT
        assert_eq!(kmers.len(), 3);
        assert_eq!(kmers[0].0, 0); // position of AC
        assert_eq!(kmers[1].0, 1); // position of CG
        assert_eq!(kmers[2].0, 2); // position of GT
    }

    #[test]
    fn test_kmer_extraction_k1() {
        let seq = encode_sequence(b"ACGT");
        let kmers: Vec<_> = KmerIterator::new(&seq, 1).collect();
        assert_eq!(kmers.len(), 4);
    }

    #[test]
    fn test_kmer_extraction_full_length() {
        let seq = encode_sequence(b"ACGT");
        let kmers: Vec<_> = KmerIterator::new(&seq, 4).collect();
        assert_eq!(kmers.len(), 1);
    }

    #[test]
    fn test_kmer_skip_n() {
        let seq = encode_sequence(b"ACNGT");
        let kmers: Vec<_> = KmerIterator::new(&seq, 2).collect();
        // AC, then N resets, then GT
        assert_eq!(kmers.len(), 2);
    }

    #[test]
    fn test_reverse_complement_kmer() {
        // AT (00 11) -> reverse complement is AT (00 11)
        let at = 0b0011u64;
        assert_eq!(reverse_complement_kmer(at, 2), at);

        // AC (00 01) -> reverse complement is GT (10 11)
        let ac = 0b0001u64;
        let gt = 0b1011u64;
        assert_eq!(reverse_complement_kmer(ac, 2), gt);
    }

    #[test]
    fn test_canonical_kmer() {
        // AC (00 01) vs GT (10 11) - AC < GT, so canonical is AC
        let ac = 0b0001u64;
        let gt = 0b1011u64;
        assert_eq!(canonical_kmer(ac, 2), ac);
        assert_eq!(canonical_kmer(gt, 2), ac);
    }

    #[test]
    fn test_kmer_to_string() {
        // AC = 00 01
        let ac = 0b0001u64;
        assert_eq!(kmer_to_string(ac, 2), "AC");

        // ACGT = 00 01 10 11
        let acgt = 0b00011011u64;
        assert_eq!(kmer_to_string(acgt, 4), "ACGT");
    }

    #[test]
    fn test_kmer_values() {
        // Verify the kmer values are what we expect
        let seq = encode_sequence(b"ACGT");
        let kmers: Vec<_> = KmerIterator::new(&seq, 2).collect();

        // AC = 00 01 = 1
        assert_eq!(kmers[0].1, 0b0001);
        // CG = 01 10 = 6
        assert_eq!(kmers[1].1, 0b0110);
        // GT = 10 11 = 11
        assert_eq!(kmers[2].1, 0b1011);
    }

    #[test]
    fn test_empty_sequence() {
        let seq: Vec<u8> = vec![];
        let kmers: Vec<_> = KmerIterator::new(&seq, 2).collect();
        assert_eq!(kmers.len(), 0);
    }

    #[test]
    fn test_sequence_shorter_than_k() {
        let seq = encode_sequence(b"AC");
        let kmers: Vec<_> = KmerIterator::new(&seq, 3).collect();
        assert_eq!(kmers.len(), 0);
    }
}
