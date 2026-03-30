pub const A: u8 = 0;
pub const C: u8 = 1;
pub const G: u8 = 2;
pub const T: u8 = 3;
pub const N: u8 = 4;

pub fn encode_base(b: u8) -> u8 {
    match b {
        b'A' | b'a' => A,
        b'C' | b'c' => C,
        b'G' | b'g' => G,
        b'T' | b't' => T,
        _ => N,
    }
}

pub fn decode_base(e: u8) -> u8 {
    match e {
        A => b'A',
        C => b'C',
        G => b'G',
        T => b'T',
        _ => b'N',
    }
}

pub fn complement(e: u8) -> u8 {
    match e {
        A => T,
        T => A,
        C => G,
        G => C,
        _ => N,
    }
}

pub fn encode_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| encode_base(b)).collect()
}

pub fn decode_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| decode_base(b)).collect()
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_decode_roundtrip() {
        assert_eq!(decode_base(encode_base(b'A')), b'A');
        assert_eq!(decode_base(encode_base(b'C')), b'C');
        assert_eq!(decode_base(encode_base(b'G')), b'G');
        assert_eq!(decode_base(encode_base(b'T')), b'T');
    }

    #[test]
    fn test_encode_case_insensitive() {
        assert_eq!(encode_base(b'a'), encode_base(b'A'));
        assert_eq!(encode_base(b'c'), encode_base(b'C'));
        assert_eq!(encode_base(b'g'), encode_base(b'G'));
        assert_eq!(encode_base(b't'), encode_base(b'T'));
    }

    #[test]
    fn test_unknown_bases_map_to_n() {
        assert_eq!(encode_base(b'X'), N);
        assert_eq!(encode_base(b'N'), N);
        assert_eq!(encode_base(b'.'), N);
    }

    #[test]
    fn test_complement() {
        assert_eq!(complement(A), T);
        assert_eq!(complement(T), A);
        assert_eq!(complement(C), G);
        assert_eq!(complement(G), C);
        assert_eq!(complement(N), N);
    }

    #[test]
    fn test_encode_sequence() {
        let seq = encode_sequence(b"ACGT");
        assert_eq!(seq, vec![A, C, G, T]);
    }

    #[test]
    fn test_reverse_complement() {
        let seq = encode_sequence(b"ACGT");
        let rc = reverse_complement(&seq);
        // ACGT -> reverse: TGCA -> complement: ACGT
        assert_eq!(rc, vec![A, C, G, T]);
    }
}
