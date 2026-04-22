/// Alignment scoring and quality assessment
///
/// This module provides tools for evaluating alignment quality,
/// computing mapping quality scores, and filtering alignments.

use super::{Cigar, CigarOp};

/// Configuration for alignment quality assessment
#[derive(Debug, Clone)]
pub struct AlignmentScorer {
    /// Minimum alignment score to consider valid (default: 30)
    pub min_score: f64,
    /// Minimum percent identity required (default: 0.80 = 80%)
    pub min_identity: f64,
    /// Maximum gap rate allowed (default: 0.10 = 10%)
    pub max_gap_rate: f64,
}

impl Default for AlignmentScorer {
    fn default() -> Self {
        Self {
            min_score: 30.0,
            min_identity: 0.80,
            max_gap_rate: 0.10,
        }
    }
}

/// Statistics computed from a CIGAR string
#[derive(Debug, Clone, Default)]
pub struct AlignmentStats {
    /// Number of matching/mismatching positions
    pub matches: u32,
    /// Number of insertions (bases in read not in reference)
    pub insertions: u32,
    /// Number of deletions (bases in reference not in read)
    pub deletions: u32,
    /// Number of soft-clipped bases
    pub soft_clips: u32,
    /// Number of hard-clipped bases
    pub hard_clips: u32,
    /// Total alignment length (M + I + D)
    pub alignment_length: u32,
}

impl AlignmentStats {
    /// Compute percent identity (matches / alignment_length)
    /// Note: Without actual sequence data, we treat all M operations as matches.
    /// In practice, you'd need the sequences to distinguish matches from mismatches.
    pub fn identity(&self) -> f64 {
        if self.alignment_length == 0 {
            return 0.0;
        }
        self.matches as f64 / self.alignment_length as f64
    }

    /// Compute gap rate ((insertions + deletions) / alignment_length)
    pub fn gap_rate(&self) -> f64 {
        if self.alignment_length == 0 {
            return 0.0;
        }
        (self.insertions + self.deletions) as f64 / self.alignment_length as f64
    }
}

/// Classification of alignment type
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentClass {
    /// Primary alignment (best hit for this read)
    Primary,
    /// Secondary alignment (alternative mapping location)
    Secondary,
    /// Supplementary alignment (chimeric alignment, part of a split read)
    Supplementary,
}

/// Result of scoring an alignment
#[derive(Debug, Clone)]
pub struct ScoredAlignment {
    /// Original alignment score
    pub score: i32,
    /// Percent identity
    pub identity: f64,
    /// Gap rate
    pub gap_rate: f64,
    /// Whether alignment passes quality thresholds
    pub passes_filter: bool,
    /// Alignment classification
    pub class: AlignmentClass,
    /// MAPQ score
    pub mapq: u8,
}

impl AlignmentScorer {
    /// Create a new scorer with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a scorer with custom thresholds
    pub fn with_thresholds(min_score: f64, min_identity: f64, max_gap_rate: f64) -> Self {
        Self {
            min_score,
            min_identity,
            max_gap_rate,
        }
    }

    /// Compute alignment statistics from a CIGAR string
    pub fn score_alignment(&self, cigar: &Cigar) -> AlignmentStats {
        let mut stats = AlignmentStats::default();

        for op in cigar.ops() {
            match op {
                CigarOp::Match(n) => {
                    stats.matches += n;
                    stats.alignment_length += n;
                }
                CigarOp::Insertion(n) => {
                    stats.insertions += n;
                    stats.alignment_length += n;
                }
                CigarOp::Deletion(n) => {
                    stats.deletions += n;
                    stats.alignment_length += n;
                }
                CigarOp::SoftClip(n) => {
                    stats.soft_clips += n;
                }
                CigarOp::HardClip(n) => {
                    stats.hard_clips += n;
                }
            }
        }

        stats
    }

    /// Check if an alignment passes quality thresholds
    pub fn is_good_alignment(&self, score: f64, identity: f64, gap_rate: f64) -> bool {
        score >= self.min_score && identity >= self.min_identity && gap_rate <= self.max_gap_rate
    }

    /// Filter alignments, keeping only those that pass quality thresholds
    pub fn filter_alignments<T: HasAlignmentScore>(&self, alignments: Vec<T>) -> Vec<T> {
        alignments
            .into_iter()
            .filter(|aln| {
                let (score, identity, gap_rate) = aln.alignment_metrics();
                self.is_good_alignment(score, identity, gap_rate)
            })
            .collect()
    }

    /// Pick the best alignment from a list
    /// Returns None if the list is empty
    /// Breaks ties by identity (higher is better)
    pub fn pick_best<T: HasAlignmentScore + Clone>(&self, alignments: &[T]) -> Option<T> {
        if alignments.is_empty() {
            return None;
        }

        let mut best_idx = 0;
        let mut best_score = f64::NEG_INFINITY;
        let mut best_identity = 0.0;

        for (i, aln) in alignments.iter().enumerate() {
            let (score, identity, _) = aln.alignment_metrics();
            if score > best_score || (score == best_score && identity > best_identity) {
                best_idx = i;
                best_score = score;
                best_identity = identity;
            }
        }

        Some(alignments[best_idx].clone())
    }

    /// Compute MAPQ from best and second-best scores
    ///
    /// MAPQ = min(60, -10 * log10(1 - best_score / (best_score + subopt_score)))
    pub fn compute_mapq(&self, best_score: i32, second_best_score: i32) -> u8 {
        if best_score <= 0 {
            return 0;
        }

        let total = best_score + second_best_score;
        if total == 0 {
            return 60;
        }

        let ratio = best_score as f64 / total as f64;
        let prob_wrong = 1.0 - ratio;

        if prob_wrong <= 0.000001 {
            return 60;
        }

        let mapq = (-10.0 * prob_wrong.log10()).min(60.0).max(0.0) as u8;
        mapq
    }

    /// Classify alignment based on CIGAR pattern
    ///
    /// - Primary: Standard alignment without major clipping
    /// - Secondary: Alternative alignment (determined by flag, not CIGAR)
    /// - Supplementary: Chimeric alignment with significant clipping (>50% clipped)
    pub fn classify_alignment(&self, cigar: &Cigar, is_secondary: bool) -> AlignmentClass {
        if is_secondary {
            return AlignmentClass::Secondary;
        }

        let stats = self.score_alignment(cigar);
        let total_length = stats.alignment_length + stats.soft_clips + stats.hard_clips;

        if total_length == 0 {
            return AlignmentClass::Primary;
        }

        let clip_rate = (stats.soft_clips + stats.hard_clips) as f64 / total_length as f64;

        // If more than 50% of the read is clipped, it's likely supplementary
        if clip_rate > 0.5 {
            AlignmentClass::Supplementary
        } else {
            AlignmentClass::Primary
        }
    }

    /// Score and evaluate an alignment comprehensively
    pub fn evaluate(
        &self,
        cigar: &Cigar,
        alignment_score: i32,
        second_best_score: i32,
        is_secondary: bool,
    ) -> ScoredAlignment {
        let stats = self.score_alignment(cigar);
        let identity = stats.identity();
        let gap_rate = stats.gap_rate();

        ScoredAlignment {
            score: alignment_score,
            identity,
            gap_rate,
            passes_filter: self.is_good_alignment(alignment_score as f64, identity, gap_rate),
            class: self.classify_alignment(cigar, is_secondary),
            mapq: self.compute_mapq(alignment_score, second_best_score),
        }
    }
}

/// Trait for types that have alignment metrics
pub trait HasAlignmentScore {
    /// Return (score, identity, gap_rate)
    fn alignment_metrics(&self) -> (f64, f64, f64);
}

/// Simple wrapper for testing
#[derive(Clone, Debug)]
pub struct SimpleAlignment {
    pub score: f64,
    pub identity: f64,
    pub gap_rate: f64,
}

impl HasAlignmentScore for SimpleAlignment {
    fn alignment_metrics(&self) -> (f64, f64, f64) {
        (self.score, self.identity, self.gap_rate)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scorer_defaults() {
        let scorer = AlignmentScorer::default();
        assert_eq!(scorer.min_score, 30.0);
        assert_eq!(scorer.min_identity, 0.80);
        assert_eq!(scorer.max_gap_rate, 0.10);
    }

    #[test]
    fn test_score_alignment_simple_match() {
        let scorer = AlignmentScorer::new();
        let cigar = Cigar::from("50M");

        let stats = scorer.score_alignment(&cigar);

        assert_eq!(stats.matches, 50);
        assert_eq!(stats.insertions, 0);
        assert_eq!(stats.deletions, 0);
        assert_eq!(stats.alignment_length, 50);
        assert_eq!(stats.identity(), 1.0);
        assert_eq!(stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_score_alignment_with_indels() {
        let scorer = AlignmentScorer::new();
        let cigar = Cigar::from("20M5I25M");

        let stats = scorer.score_alignment(&cigar);

        assert_eq!(stats.matches, 45);
        assert_eq!(stats.insertions, 5);
        assert_eq!(stats.alignment_length, 50);
        assert_eq!(stats.identity(), 0.9); // 45/50
        assert_eq!(stats.gap_rate(), 0.1); // 5/50
    }

    #[test]
    fn test_score_alignment_with_deletion() {
        let scorer = AlignmentScorer::new();
        let cigar = Cigar::from("30M10D20M");

        let stats = scorer.score_alignment(&cigar);

        assert_eq!(stats.matches, 50);
        assert_eq!(stats.deletions, 10);
        assert_eq!(stats.alignment_length, 60);
    }

    #[test]
    fn test_score_alignment_with_soft_clips() {
        let scorer = AlignmentScorer::new();
        let cigar = Cigar::from("5S40M5S");

        let stats = scorer.score_alignment(&cigar);

        assert_eq!(stats.matches, 40);
        assert_eq!(stats.soft_clips, 10);
        assert_eq!(stats.alignment_length, 40); // soft clips don't count
    }

    #[test]
    fn test_is_good_alignment() {
        let scorer = AlignmentScorer::default();

        // Good alignment: high score, high identity, low gap rate
        assert!(scorer.is_good_alignment(50.0, 0.95, 0.02));

        // Bad: low score
        assert!(!scorer.is_good_alignment(20.0, 0.95, 0.02));

        // Bad: low identity
        assert!(!scorer.is_good_alignment(50.0, 0.70, 0.02));

        // Bad: high gap rate
        assert!(!scorer.is_good_alignment(50.0, 0.95, 0.20));
    }

    #[test]
    fn test_filter_alignments() {
        let scorer = AlignmentScorer::with_thresholds(30.0, 0.80, 0.10);

        let alignments = vec![
            SimpleAlignment { score: 50.0, identity: 0.95, gap_rate: 0.02 },
            SimpleAlignment { score: 20.0, identity: 0.95, gap_rate: 0.02 }, // low score
            SimpleAlignment { score: 50.0, identity: 0.70, gap_rate: 0.02 }, // low identity
            SimpleAlignment { score: 40.0, identity: 0.85, gap_rate: 0.05 },
        ];

        let filtered = scorer.filter_alignments(alignments);

        assert_eq!(filtered.len(), 2);
        assert_eq!(filtered[0].score, 50.0);
        assert_eq!(filtered[1].score, 40.0);
    }

    #[test]
    fn test_pick_best() {
        let scorer = AlignmentScorer::new();

        let alignments = vec![
            SimpleAlignment { score: 40.0, identity: 0.90, gap_rate: 0.05 },
            SimpleAlignment { score: 50.0, identity: 0.85, gap_rate: 0.10 },
            SimpleAlignment { score: 45.0, identity: 0.95, gap_rate: 0.02 },
        ];

        let best = scorer.pick_best(&alignments).unwrap();
        assert_eq!(best.score, 50.0); // Highest score wins
    }

    #[test]
    fn test_pick_best_tie_by_identity() {
        let scorer = AlignmentScorer::new();

        let alignments = vec![
            SimpleAlignment { score: 50.0, identity: 0.85, gap_rate: 0.10 },
            SimpleAlignment { score: 50.0, identity: 0.95, gap_rate: 0.02 },
        ];

        let best = scorer.pick_best(&alignments).unwrap();
        assert_eq!(best.identity, 0.95); // Same score, higher identity wins
    }

    #[test]
    fn test_pick_best_empty() {
        let scorer = AlignmentScorer::new();
        let alignments: Vec<SimpleAlignment> = vec![];

        assert!(scorer.pick_best(&alignments).is_none());
    }

    #[test]
    fn test_compute_mapq_unique() {
        let scorer = AlignmentScorer::new();

        // Unique alignment (no second best)
        let mapq = scorer.compute_mapq(100, 0);
        assert_eq!(mapq, 60); // Maximum MAPQ for unique alignment
    }

    #[test]
    fn test_compute_mapq_with_competitor() {
        let scorer = AlignmentScorer::new();

        // Close second best
        let mapq = scorer.compute_mapq(100, 90);
        assert!(mapq < 60);
        assert!(mapq > 0);
    }

    #[test]
    fn test_compute_mapq_zero_score() {
        let scorer = AlignmentScorer::new();
        assert_eq!(scorer.compute_mapq(0, 0), 0);
        assert_eq!(scorer.compute_mapq(-10, 0), 0);
    }

    #[test]
    fn test_compute_mapq_equal_scores() {
        let scorer = AlignmentScorer::new();
        let mapq = scorer.compute_mapq(50, 50);
        // With equal scores, probability of wrong mapping is 0.5
        // MAPQ = -10 * log10(0.5) = 3.01...
        assert!(mapq >= 3 && mapq <= 4);
    }

    #[test]
    fn test_classify_alignment_primary() {
        let scorer = AlignmentScorer::new();
        let cigar = Cigar::from("100M");

        assert_eq!(
            scorer.classify_alignment(&cigar, false),
            AlignmentClass::Primary
        );
    }

    #[test]
    fn test_classify_alignment_secondary() {
        let scorer = AlignmentScorer::new();
        let cigar = Cigar::from("100M");

        assert_eq!(
            scorer.classify_alignment(&cigar, true),
            AlignmentClass::Secondary
        );
    }

    #[test]
    fn test_classify_alignment_supplementary() {
        let scorer = AlignmentScorer::new();
        // More than 50% clipped -> supplementary
        let cigar = Cigar::from("60S40M");

        assert_eq!(
            scorer.classify_alignment(&cigar, false),
            AlignmentClass::Supplementary
        );
    }

    #[test]
    fn test_classify_alignment_with_small_clips() {
        let scorer = AlignmentScorer::new();
        // Less than 50% clipped -> primary
        let cigar = Cigar::from("5S90M5S");

        assert_eq!(
            scorer.classify_alignment(&cigar, false),
            AlignmentClass::Primary
        );
    }

    #[test]
    fn test_evaluate_full() {
        let scorer = AlignmentScorer::default();
        let cigar = Cigar::from("45M5I");

        let result = scorer.evaluate(&cigar, 80, 20, false);

        assert_eq!(result.score, 80);
        assert_eq!(result.identity, 0.9); // 45/50
        assert_eq!(result.gap_rate, 0.1); // 5/50
        assert!(result.passes_filter);
        assert_eq!(result.class, AlignmentClass::Primary);
        assert!(result.mapq > 0);
    }

    #[test]
    fn test_stats_identity_empty() {
        let stats = AlignmentStats::default();
        assert_eq!(stats.identity(), 0.0);
        assert_eq!(stats.gap_rate(), 0.0);
    }

    #[test]
    fn test_scorer_with_thresholds() {
        let scorer = AlignmentScorer::with_thresholds(50.0, 0.90, 0.05);

        assert_eq!(scorer.min_score, 50.0);
        assert_eq!(scorer.min_identity, 0.90);
        assert_eq!(scorer.max_gap_rate, 0.05);

        // Should be stricter
        assert!(!scorer.is_good_alignment(40.0, 0.95, 0.02)); // score too low
        assert!(!scorer.is_good_alignment(60.0, 0.85, 0.02)); // identity too low
        assert!(!scorer.is_good_alignment(60.0, 0.95, 0.08)); // gap rate too high
        assert!(scorer.is_good_alignment(60.0, 0.95, 0.02)); // all good
    }
}
