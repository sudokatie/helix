use crate::index::{HashIndex, MinimizerIterator};
use crate::seq::Reference;

/// Configuration for index building
#[derive(Debug, Clone)]
pub struct IndexConfig {
    /// K-mer size (default: 15)
    pub k: usize,
    /// Minimizer window size (default: 10)
    pub w: usize,
}

impl Default for IndexConfig {
    fn default() -> Self {
        Self { k: 15, w: 10 }
    }
}

impl IndexConfig {
    pub fn new(k: usize, w: usize) -> Self {
        assert!(k > 0 && k <= 32, "k must be between 1 and 32");
        assert!(w > 0, "w must be positive");
        Self { k, w }
    }
}

/// Information about an indexed reference sequence
#[derive(Debug, Clone)]
pub struct ReferenceInfo {
    pub name: String,
    pub length: usize,
    /// Global position offset for this reference
    pub offset: u64,
}

/// Complete index for alignment
pub struct Index {
    pub references: Vec<ReferenceInfo>,
    pub hash: HashIndex,
    pub config: IndexConfig,
}

impl Index {
    /// Get the reference containing a global position
    pub fn position_to_reference(&self, pos: u64) -> Option<(&ReferenceInfo, u64)> {
        for (i, r) in self.references.iter().enumerate() {
            let end = if i + 1 < self.references.len() {
                self.references[i + 1].offset
            } else {
                r.offset + r.length as u64
            };

            if pos >= r.offset && pos < end {
                return Some((r, pos - r.offset));
            }
        }
        None
    }

    /// Total length of all references
    pub fn total_length(&self) -> u64 {
        self.references.iter().map(|r| r.length as u64).sum()
    }
}

/// Builder for constructing an index
pub struct IndexBuilder {
    config: IndexConfig,
    references: Vec<ReferenceInfo>,
    hash: HashIndex,
    current_offset: u64,
}

impl IndexBuilder {
    pub fn new(config: IndexConfig) -> Self {
        Self {
            config,
            references: Vec::new(),
            hash: HashIndex::new(),
            current_offset: 0,
        }
    }

    /// Add a reference sequence to the index
    pub fn add_reference(&mut self, reference: &Reference) {
        let info = ReferenceInfo {
            name: reference.name.clone(),
            length: reference.seq.len(),
            offset: self.current_offset,
        };

        // Extract minimizers and add to hash
        for minimizer in MinimizerIterator::new(&reference.seq, self.config.k, self.config.w) {
            let global_pos = self.current_offset + minimizer.pos as u64;
            self.hash.insert(minimizer.kmer, global_pos as u32);
        }

        self.current_offset += reference.seq.len() as u64;
        self.references.push(info);
    }

    /// Build the final index
    pub fn build(mut self) -> Index {
        self.hash.finalize();
        Index {
            references: self.references,
            hash: self.hash,
            config: self.config,
        }
    }

    /// Number of references added so far
    pub fn reference_count(&self) -> usize {
        self.references.len()
    }

    /// Total bases indexed so far
    pub fn total_bases(&self) -> u64 {
        self.current_offset
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::seq::encode_sequence;

    fn make_reference(name: &str, seq: &[u8]) -> Reference {
        Reference {
            name: name.to_string(),
            seq: encode_sequence(seq),
        }
    }

    #[test]
    fn test_build_single_reference() {
        let config = IndexConfig::new(3, 2);
        let mut builder = IndexBuilder::new(config);

        let r = make_reference("chr1", b"ACGTACGTACGT");
        builder.add_reference(&r);

        let index = builder.build();

        assert_eq!(index.references.len(), 1);
        assert_eq!(index.references[0].name, "chr1");
        assert_eq!(index.references[0].length, 12);
        assert_eq!(index.references[0].offset, 0);
    }

    #[test]
    fn test_build_multiple_references() {
        let config = IndexConfig::new(3, 2);
        let mut builder = IndexBuilder::new(config);

        builder.add_reference(&make_reference("chr1", b"ACGTACGT"));
        builder.add_reference(&make_reference("chr2", b"TGCATGCA"));

        let index = builder.build();

        assert_eq!(index.references.len(), 2);
        assert_eq!(index.references[0].offset, 0);
        assert_eq!(index.references[1].offset, 8);
        assert_eq!(index.total_length(), 16);
    }

    #[test]
    fn test_position_to_reference() {
        let config = IndexConfig::new(3, 2);
        let mut builder = IndexBuilder::new(config);

        builder.add_reference(&make_reference("chr1", b"ACGTACGT"));
        builder.add_reference(&make_reference("chr2", b"TGCATGCA"));

        let index = builder.build();

        // Position in chr1
        let (r, local) = index.position_to_reference(5).unwrap();
        assert_eq!(r.name, "chr1");
        assert_eq!(local, 5);

        // Position in chr2
        let (r, local) = index.position_to_reference(10).unwrap();
        assert_eq!(r.name, "chr2");
        assert_eq!(local, 2);

        // Position beyond index
        assert!(index.position_to_reference(100).is_none());
    }

    #[test]
    fn test_minimizers_indexed() {
        let config = IndexConfig::new(3, 2);
        let mut builder = IndexBuilder::new(config);

        builder.add_reference(&make_reference("chr1", b"ACGTACGTACGT"));
        let index = builder.build();

        // Index should have some minimizers
        assert!(!index.hash.is_empty());
    }

    #[test]
    fn test_config_defaults() {
        let config = IndexConfig::default();
        assert_eq!(config.k, 15);
        assert_eq!(config.w, 10);
    }

    #[test]
    fn test_builder_stats() {
        let config = IndexConfig::new(3, 2);
        let mut builder = IndexBuilder::new(config);

        assert_eq!(builder.reference_count(), 0);
        assert_eq!(builder.total_bases(), 0);

        builder.add_reference(&make_reference("chr1", b"ACGTACGT"));
        assert_eq!(builder.reference_count(), 1);
        assert_eq!(builder.total_bases(), 8);
    }
}
