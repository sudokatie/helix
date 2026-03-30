use crate::index::kmer::Kmer;
use std::collections::HashMap;

/// Hash index for k-mer to position lookups
pub struct HashIndex {
    // During building: kmer -> positions
    builder: Option<HashMap<Kmer, Vec<u32>>>,
    // After finalize: flattened for cache efficiency
    table: Vec<u64>,         // kmer values (sorted)
    offsets: Vec<u32>,       // Start offset for each kmer's positions
    positions: Vec<u32>,     // All positions, grouped by kmer
}

impl HashIndex {
    pub fn new() -> Self {
        Self {
            builder: Some(HashMap::new()),
            table: Vec::new(),
            offsets: Vec::new(),
            positions: Vec::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            builder: Some(HashMap::with_capacity(capacity)),
            table: Vec::new(),
            offsets: Vec::new(),
            positions: Vec::new(),
        }
    }

    /// Insert a k-mer at a position (during building phase)
    pub fn insert(&mut self, kmer: Kmer, pos: u32) {
        if let Some(ref mut builder) = self.builder {
            builder.entry(kmer).or_insert_with(Vec::new).push(pos);
        } else {
            panic!("Cannot insert after finalize");
        }
    }

    /// Finalize the index for lookups (sorts and flattens)
    pub fn finalize(&mut self) {
        let builder = self.builder.take().expect("Already finalized");

        // Sort k-mers for binary search
        let mut kmers: Vec<_> = builder.into_iter().collect();
        kmers.sort_by_key(|(k, _)| *k);

        self.table = Vec::with_capacity(kmers.len());
        self.offsets = Vec::with_capacity(kmers.len() + 1);
        
        let total_positions: usize = kmers.iter().map(|(_, v)| v.len()).sum();
        self.positions = Vec::with_capacity(total_positions);

        for (kmer, mut positions) in kmers {
            self.table.push(kmer);
            self.offsets.push(self.positions.len() as u32);
            positions.sort();
            self.positions.extend(positions);
        }
        self.offsets.push(self.positions.len() as u32);
    }

    /// Lookup positions for a k-mer
    pub fn lookup(&self, kmer: Kmer) -> &[u32] {
        if self.builder.is_some() {
            panic!("Must call finalize() before lookup");
        }

        match self.table.binary_search(&kmer) {
            Ok(idx) => {
                let start = self.offsets[idx] as usize;
                let end = self.offsets[idx + 1] as usize;
                &self.positions[start..end]
            }
            Err(_) => &[],
        }
    }

    /// Number of unique k-mers in the index
    pub fn len(&self) -> usize {
        if let Some(ref builder) = self.builder {
            builder.len()
        } else {
            self.table.len()
        }
    }

    /// Check if index is empty
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Total number of positions stored
    pub fn total_positions(&self) -> usize {
        if let Some(ref builder) = self.builder {
            builder.values().map(|v| v.len()).sum()
        } else {
            self.positions.len()
        }
    }

    /// Check if index has been finalized
    pub fn is_finalized(&self) -> bool {
        self.builder.is_none()
    }
}

impl Default for HashIndex {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert_and_lookup_single() {
        let mut index = HashIndex::new();
        index.insert(42, 100);
        index.finalize();

        let positions = index.lookup(42);
        assert_eq!(positions, &[100]);
    }

    #[test]
    fn test_multiple_positions_same_kmer() {
        let mut index = HashIndex::new();
        index.insert(42, 100);
        index.insert(42, 200);
        index.insert(42, 150);
        index.finalize();

        let positions = index.lookup(42);
        assert_eq!(positions, &[100, 150, 200]); // Sorted
    }

    #[test]
    fn test_lookup_missing_kmer() {
        let mut index = HashIndex::new();
        index.insert(42, 100);
        index.finalize();

        let positions = index.lookup(99);
        assert!(positions.is_empty());
    }

    #[test]
    fn test_multiple_kmers() {
        let mut index = HashIndex::new();
        index.insert(10, 1);
        index.insert(20, 2);
        index.insert(30, 3);
        index.insert(20, 4);
        index.finalize();

        assert_eq!(index.lookup(10), &[1]);
        assert_eq!(index.lookup(20), &[2, 4]);
        assert_eq!(index.lookup(30), &[3]);
        assert_eq!(index.lookup(40), &[]);
    }

    #[test]
    fn test_large_index() {
        let mut index = HashIndex::with_capacity(10000);
        
        // Insert many k-mers
        for i in 0..10000u64 {
            index.insert(i, (i * 100) as u32);
            if i % 3 == 0 {
                index.insert(i, (i * 100 + 1) as u32);
            }
        }
        index.finalize();

        // Verify some lookups
        assert_eq!(index.lookup(0), &[0, 1]);
        assert_eq!(index.lookup(1), &[100]);
        assert_eq!(index.lookup(3), &[300, 301]);
        assert_eq!(index.lookup(9999), &[999900, 999901]);
    }

    #[test]
    fn test_len_and_total_positions() {
        let mut index = HashIndex::new();
        index.insert(10, 1);
        index.insert(20, 2);
        index.insert(20, 3);
        
        assert_eq!(index.len(), 2);
        assert_eq!(index.total_positions(), 3);
        
        index.finalize();
        
        assert_eq!(index.len(), 2);
        assert_eq!(index.total_positions(), 3);
    }

    #[test]
    fn test_empty_index() {
        let mut index = HashIndex::new();
        assert!(index.is_empty());
        index.finalize();
        assert!(index.is_empty());
        assert_eq!(index.lookup(42), &[]);
    }

    #[test]
    fn test_is_finalized() {
        let mut index = HashIndex::new();
        assert!(!index.is_finalized());
        index.finalize();
        assert!(index.is_finalized());
    }
}
