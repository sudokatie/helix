use std::fmt;

/// CIGAR operation types
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    /// Alignment match (can be match or mismatch)
    Match(u32),
    /// Insertion to the reference (gap in reference)
    Insertion(u32),
    /// Deletion from the reference (gap in query)
    Deletion(u32),
    /// Soft clipping (unaligned at ends, sequence present)
    SoftClip(u32),
    /// Hard clipping (unaligned at ends, sequence not present)
    HardClip(u32),
}

impl CigarOp {
    /// Get the operation character
    pub fn op_char(&self) -> char {
        match self {
            CigarOp::Match(_) => 'M',
            CigarOp::Insertion(_) => 'I',
            CigarOp::Deletion(_) => 'D',
            CigarOp::SoftClip(_) => 'S',
            CigarOp::HardClip(_) => 'H',
        }
    }

    /// Get the length of this operation
    pub fn len(&self) -> u32 {
        match self {
            CigarOp::Match(n)
            | CigarOp::Insertion(n)
            | CigarOp::Deletion(n)
            | CigarOp::SoftClip(n)
            | CigarOp::HardClip(n) => *n,
        }
    }

    /// Check if operation is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Does this operation consume query bases?
    pub fn consumes_query(&self) -> bool {
        matches!(
            self,
            CigarOp::Match(_) | CigarOp::Insertion(_) | CigarOp::SoftClip(_)
        )
    }

    /// Does this operation consume reference bases?
    pub fn consumes_reference(&self) -> bool {
        matches!(self, CigarOp::Match(_) | CigarOp::Deletion(_))
    }
}

impl fmt::Display for CigarOp {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.len(), self.op_char())
    }
}

/// CIGAR string representing an alignment
#[derive(Clone, Debug, Default)]
pub struct Cigar {
    ops: Vec<CigarOp>,
}

impl Cigar {
    /// Create a new empty CIGAR
    pub fn new() -> Self {
        Self { ops: Vec::new() }
    }

    /// Add an operation, merging with previous if same type
    pub fn push(&mut self, op: CigarOp) {
        if op.is_empty() {
            return;
        }

        if let Some(last) = self.ops.last_mut() {
            // Merge consecutive operations of the same type
            let merged = match (last, &op) {
                (CigarOp::Match(n), CigarOp::Match(m)) => {
                    *n += m;
                    true
                }
                (CigarOp::Insertion(n), CigarOp::Insertion(m)) => {
                    *n += m;
                    true
                }
                (CigarOp::Deletion(n), CigarOp::Deletion(m)) => {
                    *n += m;
                    true
                }
                (CigarOp::SoftClip(n), CigarOp::SoftClip(m)) => {
                    *n += m;
                    true
                }
                (CigarOp::HardClip(n), CigarOp::HardClip(m)) => {
                    *n += m;
                    true
                }
                _ => false,
            };
            if merged {
                return;
            }
        }

        self.ops.push(op);
    }

    /// Add a single match/mismatch
    pub fn push_match(&mut self) {
        self.push(CigarOp::Match(1));
    }

    /// Add a single insertion
    pub fn push_insertion(&mut self) {
        self.push(CigarOp::Insertion(1));
    }

    /// Add a single deletion
    pub fn push_deletion(&mut self) {
        self.push(CigarOp::Deletion(1));
    }

    /// Get the operations
    pub fn ops(&self) -> &[CigarOp] {
        &self.ops
    }

    /// Check if empty
    pub fn is_empty(&self) -> bool {
        self.ops.is_empty()
    }

    /// Number of operations
    pub fn len(&self) -> usize {
        self.ops.len()
    }

    /// Total query bases consumed
    pub fn query_length(&self) -> u32 {
        self.ops
            .iter()
            .filter(|op| op.consumes_query())
            .map(|op| op.len())
            .sum()
    }

    /// Total reference bases consumed
    pub fn reference_length(&self) -> u32 {
        self.ops
            .iter()
            .filter(|op| op.consumes_reference())
            .map(|op| op.len())
            .sum()
    }

    /// Reverse the CIGAR (for traceback)
    pub fn reverse(&mut self) {
        self.ops.reverse();
    }
}

impl fmt::Display for Cigar {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.ops.is_empty() {
            return write!(f, "*");
        }
        for op in &self.ops {
            write!(f, "{}", op)?;
        }
        Ok(())
    }
}

impl From<&str> for Cigar {
    fn from(s: &str) -> Self {
        let mut cigar = Cigar::new();
        let mut num = 0u32;

        for c in s.chars() {
            if c.is_ascii_digit() {
                num = num * 10 + c.to_digit(10).unwrap();
            } else {
                let op = match c {
                    'M' => CigarOp::Match(num),
                    'I' => CigarOp::Insertion(num),
                    'D' => CigarOp::Deletion(num),
                    'S' => CigarOp::SoftClip(num),
                    'H' => CigarOp::HardClip(num),
                    _ => continue,
                };
                if num > 0 {
                    cigar.ops.push(op);
                }
                num = 0;
            }
        }

        cigar
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cigar_match() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Match(10));
        assert_eq!(cigar.to_string(), "10M");
    }

    #[test]
    fn test_cigar_merge_consecutive() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Match(5));
        cigar.push(CigarOp::Match(5));
        assert_eq!(cigar.to_string(), "10M");
        assert_eq!(cigar.len(), 1);
    }

    #[test]
    fn test_cigar_with_insertion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Match(5));
        cigar.push(CigarOp::Insertion(2));
        cigar.push(CigarOp::Match(3));
        assert_eq!(cigar.to_string(), "5M2I3M");
    }

    #[test]
    fn test_cigar_with_deletion() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::Match(5));
        cigar.push(CigarOp::Deletion(2));
        cigar.push(CigarOp::Match(3));
        assert_eq!(cigar.to_string(), "5M2D3M");
    }

    #[test]
    fn test_cigar_soft_clips() {
        let mut cigar = Cigar::new();
        cigar.push(CigarOp::SoftClip(3));
        cigar.push(CigarOp::Match(10));
        cigar.push(CigarOp::SoftClip(2));
        assert_eq!(cigar.to_string(), "3S10M2S");
    }

    #[test]
    fn test_cigar_parse() {
        let cigar = Cigar::from("5M2I3M");
        assert_eq!(cigar.ops().len(), 3);
        assert_eq!(cigar.ops()[0], CigarOp::Match(5));
        assert_eq!(cigar.ops()[1], CigarOp::Insertion(2));
        assert_eq!(cigar.ops()[2], CigarOp::Match(3));
    }

    #[test]
    fn test_cigar_query_length() {
        let cigar = Cigar::from("5M2I3M2D4M");
        // M consumes query: 5 + 3 + 4 = 12
        // I consumes query: 2
        // D does not consume query
        assert_eq!(cigar.query_length(), 14);
    }

    #[test]
    fn test_cigar_reference_length() {
        let cigar = Cigar::from("5M2I3M2D4M");
        // M consumes reference: 5 + 3 + 4 = 12
        // D consumes reference: 2
        // I does not consume reference
        assert_eq!(cigar.reference_length(), 14);
    }

    #[test]
    fn test_empty_cigar() {
        let cigar = Cigar::new();
        assert_eq!(cigar.to_string(), "*");
        assert!(cigar.is_empty());
    }

    #[test]
    fn test_cigar_reverse() {
        let mut cigar = Cigar::from("5M2I3M");
        cigar.reverse();
        assert_eq!(cigar.to_string(), "3M2I5M");
    }

    #[test]
    fn test_push_helpers() {
        let mut cigar = Cigar::new();
        cigar.push_match();
        cigar.push_match();
        cigar.push_insertion();
        cigar.push_match();
        assert_eq!(cigar.to_string(), "2M1I1M");
    }
}
