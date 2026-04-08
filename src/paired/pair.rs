// Paired-end read alignment
// Handles concordant pair detection and proper SAM flag setting

use crate::output::{
    SamRecord, FLAG_MATE_REVERSE, FLAG_MATE_UNMAPPED, FLAG_PAIRED, FLAG_PROPER_PAIR,
    FLAG_READ1, FLAG_READ2, FLAG_REVERSE, FLAG_UNMAPPED,
};
use crate::paired::insert::{InsertSizeEstimator, InsertSizeStats};
use crate::seq::Sequence;

/// Orientation of a paired-end read pair
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PairOrientation {
    /// Forward-Reverse (most common for Illumina)
    FR,
    /// Reverse-Forward
    RF,
    /// Forward-Forward
    FF,
    /// Reverse-Reverse
    RR,
}

impl PairOrientation {
    /// Detect orientation from alignment strands
    pub fn from_strands(read1_reverse: bool, read2_reverse: bool) -> Self {
        match (read1_reverse, read2_reverse) {
            (false, true) => PairOrientation::FR,
            (true, false) => PairOrientation::RF,
            (false, false) => PairOrientation::FF,
            (true, true) => PairOrientation::RR,
        }
    }

    /// Check if this is the expected orientation for standard Illumina libraries
    pub fn is_expected_illumina(&self) -> bool {
        matches!(self, PairOrientation::FR | PairOrientation::RF)
    }
}

/// A pair of reads
#[derive(Debug, Clone)]
pub struct ReadPair {
    pub read1: Sequence,
    pub read2: Sequence,
}

impl ReadPair {
    pub fn new(read1: Sequence, read2: Sequence) -> Self {
        Self { read1, read2 }
    }
}

/// Alignment result for a read
#[derive(Debug, Clone)]
pub struct ReadAlignment {
    /// Reference name
    pub rname: String,
    /// 0-based position on reference
    pub pos: u32,
    /// Is reverse strand
    pub is_reverse: bool,
    /// CIGAR string
    pub cigar: String,
    /// Mapping quality
    pub mapq: u8,
    /// Is this read mapped
    pub is_mapped: bool,
}

/// Result of aligning a pair
#[derive(Debug, Clone)]
pub struct PairAlignmentResult {
    pub read1_alignment: Option<ReadAlignment>,
    pub read2_alignment: Option<ReadAlignment>,
    pub is_proper_pair: bool,
    pub insert_size: i32,
    pub orientation: Option<PairOrientation>,
}

impl PairAlignmentResult {
    /// Create an unmapped pair result
    pub fn unmapped() -> Self {
        Self {
            read1_alignment: None,
            read2_alignment: None,
            is_proper_pair: false,
            insert_size: 0,
            orientation: None,
        }
    }

    /// Calculate insert size from alignments
    /// Template length = rightmost end - leftmost start + 1
    pub fn calculate_insert_size(
        pos1: u32,
        len1: usize,
        is_rev1: bool,
        pos2: u32,
        len2: usize,
        is_rev2: bool,
    ) -> i32 {
        let end1 = pos1 as i32 + len1 as i32;
        let end2 = pos2 as i32 + len2 as i32;
        let start1 = pos1 as i32;
        let start2 = pos2 as i32;

        let leftmost = start1.min(start2);
        let rightmost = end1.max(end2);
        let tlen = rightmost - leftmost;

        // Convention: positive for read1 < read2, negative otherwise
        if (start1 <= start2 && !is_rev1) || (start2 <= start1 && is_rev2) {
            tlen
        } else {
            -tlen
        }
    }
}

/// Paired-end aligner that handles concordant pair detection
pub struct PairedAligner {
    /// Expected orientation (default FR for Illumina)
    pub expected_orientation: PairOrientation,
    /// Insert size estimator
    insert_estimator: InsertSizeEstimator,
    /// Fixed insert size stats (if provided)
    fixed_insert_stats: Option<InsertSizeStats>,
}

impl PairedAligner {
    pub fn new() -> Self {
        Self {
            expected_orientation: PairOrientation::FR,
            insert_estimator: InsertSizeEstimator::new(),
            fixed_insert_stats: None,
        }
    }

    /// Set expected orientation
    pub fn with_orientation(mut self, orientation: PairOrientation) -> Self {
        self.expected_orientation = orientation;
        self
    }

    /// Set fixed insert size statistics (skip estimation)
    pub fn with_insert_stats(mut self, stats: InsertSizeStats) -> Self {
        self.fixed_insert_stats = Some(stats);
        self
    }

    /// Get current insert size statistics
    pub fn insert_stats(&self) -> InsertSizeStats {
        self.fixed_insert_stats
            .clone()
            .unwrap_or_else(|| self.insert_estimator.stats())
    }

    /// Check if a pair is concordant
    pub fn is_concordant(
        &self,
        aln1: &ReadAlignment,
        aln2: &ReadAlignment,
        read1_len: usize,
        read2_len: usize,
    ) -> bool {
        // Must be on same reference
        if aln1.rname != aln2.rname {
            return false;
        }

        // Check orientation
        let orientation = PairOrientation::from_strands(aln1.is_reverse, aln2.is_reverse);
        if !matches!(
            (self.expected_orientation, orientation),
            (PairOrientation::FR, PairOrientation::FR)
                | (PairOrientation::FR, PairOrientation::RF)
                | (PairOrientation::RF, PairOrientation::RF)
                | (PairOrientation::RF, PairOrientation::FR)
        ) {
            // For FR library: read1 should be forward, read2 reverse (or vice versa)
            // This is a simplified check
            if !orientation.is_expected_illumina() {
                return false;
            }
        }

        // Check insert size
        let insert_size = PairAlignmentResult::calculate_insert_size(
            aln1.pos,
            read1_len,
            aln1.is_reverse,
            aln2.pos,
            read2_len,
            aln2.is_reverse,
        );

        let stats = self.insert_stats();
        stats.is_concordant(insert_size.abs())
    }

    /// Record insert size from a concordant pair
    pub fn record_insert_size(&mut self, insert_size: i32) {
        if self.fixed_insert_stats.is_none() {
            self.insert_estimator.add_sample(insert_size.abs());
        }
    }

    /// Convert alignment result to SAM records for a pair
    pub fn to_sam_records(
        &self,
        pair: &ReadPair,
        result: &PairAlignmentResult,
    ) -> (SamRecord, SamRecord) {
        let mut rec1 = if let Some(ref aln) = result.read1_alignment {
            SamRecord {
                qname: pair.read1.id.clone(),
                flag: FLAG_PAIRED | FLAG_READ1,
                rname: aln.rname.clone(),
                pos: aln.pos + 1, // SAM is 1-based
                mapq: aln.mapq,
                cigar: aln.cigar.clone(),
                rnext: "=".to_string(),
                pnext: 0,
                tlen: result.insert_size,
                seq: String::from_utf8_lossy(&decode_sequence(&pair.read1.seq)).to_string(),
                qual: String::from_utf8_lossy(&pair.read1.qual).to_string(),
                tags: Vec::new(),
            }
        } else {
            SamRecord::new_unmapped(
                pair.read1.id.clone(),
                String::from_utf8_lossy(&decode_sequence(&pair.read1.seq)).to_string(),
                String::from_utf8_lossy(&pair.read1.qual).to_string(),
            )
        };

        let mut rec2 = if let Some(ref aln) = result.read2_alignment {
            SamRecord {
                qname: pair.read2.id.clone(),
                flag: FLAG_PAIRED | FLAG_READ2,
                rname: aln.rname.clone(),
                pos: aln.pos + 1,
                mapq: aln.mapq,
                cigar: aln.cigar.clone(),
                rnext: "=".to_string(),
                pnext: 0,
                tlen: -result.insert_size, // Negative for read2
                seq: String::from_utf8_lossy(&decode_sequence(&pair.read2.seq)).to_string(),
                qual: String::from_utf8_lossy(&pair.read2.qual).to_string(),
                tags: Vec::new(),
            }
        } else {
            SamRecord::new_unmapped(
                pair.read2.id.clone(),
                String::from_utf8_lossy(&decode_sequence(&pair.read2.seq)).to_string(),
                String::from_utf8_lossy(&pair.read2.qual).to_string(),
            )
        };

        // Set mate information
        self.set_mate_info(&mut rec1, &mut rec2, result);

        (rec1, rec2)
    }

    fn set_mate_info(&self, rec1: &mut SamRecord, rec2: &mut SamRecord, result: &PairAlignmentResult) {
        let aln1 = &result.read1_alignment;
        let aln2 = &result.read2_alignment;

        // Set proper pair flag
        if result.is_proper_pair {
            rec1.flag |= FLAG_PROPER_PAIR;
            rec2.flag |= FLAG_PROPER_PAIR;
        }

        // Set strand flags
        if let Some(ref a1) = aln1 {
            if a1.is_reverse {
                rec1.flag |= FLAG_REVERSE;
            }
        }
        if let Some(ref a2) = aln2 {
            if a2.is_reverse {
                rec2.flag |= FLAG_REVERSE;
            }
        }

        // Set mate strand flags
        if let Some(ref a2) = aln2 {
            if a2.is_reverse {
                rec1.flag |= FLAG_MATE_REVERSE;
            }
            rec1.pnext = a2.pos + 1;
            rec1.rnext = if aln1.as_ref().map(|a| &a.rname) == Some(&a2.rname) {
                "=".to_string()
            } else {
                a2.rname.clone()
            };
        } else {
            rec1.flag |= FLAG_MATE_UNMAPPED;
            rec1.rnext = "*".to_string();
            rec1.pnext = 0;
        }

        if let Some(ref a1) = aln1 {
            if a1.is_reverse {
                rec2.flag |= FLAG_MATE_REVERSE;
            }
            rec2.pnext = a1.pos + 1;
            rec2.rnext = if aln2.as_ref().map(|a| &a.rname) == Some(&a1.rname) {
                "=".to_string()
            } else {
                a1.rname.clone()
            };
        } else {
            rec2.flag |= FLAG_MATE_UNMAPPED;
            rec2.rnext = "*".to_string();
            rec2.pnext = 0;
        }

        // Handle unmapped reads
        if aln1.is_none() {
            rec1.flag |= FLAG_UNMAPPED;
        }
        if aln2.is_none() {
            rec2.flag |= FLAG_UNMAPPED;
        }
    }
}

impl Default for PairedAligner {
    fn default() -> Self {
        Self::new()
    }
}

/// Decode encoded sequence back to ASCII
fn decode_sequence(encoded: &[u8]) -> Vec<u8> {
    const DECODE: [u8; 5] = [b'A', b'C', b'G', b'T', b'N'];
    encoded
        .iter()
        .map(|&b| DECODE.get(b as usize).copied().unwrap_or(b'N'))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pair_orientation() {
        assert_eq!(
            PairOrientation::from_strands(false, true),
            PairOrientation::FR
        );
        assert_eq!(
            PairOrientation::from_strands(true, false),
            PairOrientation::RF
        );
        assert_eq!(
            PairOrientation::from_strands(false, false),
            PairOrientation::FF
        );
        assert_eq!(
            PairOrientation::from_strands(true, true),
            PairOrientation::RR
        );
    }

    #[test]
    fn test_expected_orientation() {
        assert!(PairOrientation::FR.is_expected_illumina());
        assert!(PairOrientation::RF.is_expected_illumina());
        assert!(!PairOrientation::FF.is_expected_illumina());
        assert!(!PairOrientation::RR.is_expected_illumina());
    }

    #[test]
    fn test_insert_size_calculation() {
        // Read1 at pos 100, len 150, forward
        // Read2 at pos 300, len 150, reverse
        // Template = 300 + 150 - 100 = 350
        let tlen = PairAlignmentResult::calculate_insert_size(100, 150, false, 300, 150, true);
        assert_eq!(tlen, 350);
    }

    #[test]
    fn test_paired_aligner_new() {
        let aligner = PairedAligner::new();
        assert_eq!(aligner.expected_orientation, PairOrientation::FR);
    }

    #[test]
    fn test_concordant_detection() {
        let aligner = PairedAligner::new().with_insert_stats(InsertSizeStats {
            mean: 300.0,
            std_dev: 50.0,
            min: 100,
            max: 500,
            count: 1000,
            median: 300,
        });

        let aln1 = ReadAlignment {
            rname: "chr1".to_string(),
            pos: 100,
            is_reverse: false,
            cigar: "150M".to_string(),
            mapq: 60,
            is_mapped: true,
        };

        let aln2 = ReadAlignment {
            rname: "chr1".to_string(),
            pos: 300,
            is_reverse: true,
            cigar: "150M".to_string(),
            mapq: 60,
            is_mapped: true,
        };

        assert!(aligner.is_concordant(&aln1, &aln2, 150, 150));
    }

    #[test]
    fn test_discordant_different_chr() {
        let aligner = PairedAligner::new().with_insert_stats(InsertSizeStats::default());

        let aln1 = ReadAlignment {
            rname: "chr1".to_string(),
            pos: 100,
            is_reverse: false,
            cigar: "150M".to_string(),
            mapq: 60,
            is_mapped: true,
        };

        let aln2 = ReadAlignment {
            rname: "chr2".to_string(),
            pos: 100,
            is_reverse: true,
            cigar: "150M".to_string(),
            mapq: 60,
            is_mapped: true,
        };

        assert!(!aligner.is_concordant(&aln1, &aln2, 150, 150));
    }
}
