// Insert size estimation for paired-end reads
// Estimates the distribution of template lengths from initial alignments

use std::collections::VecDeque;

/// Statistics for insert size distribution
#[derive(Debug, Clone)]
pub struct InsertSizeStats {
    /// Mean insert size
    pub mean: f64,
    /// Standard deviation
    pub std_dev: f64,
    /// Minimum observed
    pub min: i32,
    /// Maximum observed
    pub max: i32,
    /// Sample count
    pub count: usize,
    /// Median insert size
    pub median: i32,
}

impl Default for InsertSizeStats {
    fn default() -> Self {
        Self {
            mean: 300.0,      // Common default for Illumina
            std_dev: 50.0,
            min: 0,
            max: 1000,
            count: 0,
            median: 300,
        }
    }
}

impl InsertSizeStats {
    /// Check if an insert size is within expected range (mean +/- 3*std_dev)
    pub fn is_expected(&self, insert_size: i32) -> bool {
        let lower = (self.mean - 3.0 * self.std_dev) as i32;
        let upper = (self.mean + 3.0 * self.std_dev) as i32;
        insert_size >= lower && insert_size <= upper
    }

    /// Check if a pair is likely concordant based on insert size
    pub fn is_concordant(&self, insert_size: i32) -> bool {
        // More lenient than is_expected, allows 6*std_dev
        let lower = (self.mean - 6.0 * self.std_dev).max(0.0) as i32;
        let upper = (self.mean + 6.0 * self.std_dev) as i32;
        insert_size >= lower && insert_size <= upper
    }
}

/// Online insert size estimator using reservoir sampling
pub struct InsertSizeEstimator {
    /// Sample buffer (limited size for memory)
    samples: VecDeque<i32>,
    /// Maximum samples to keep
    max_samples: usize,
    /// Running sum for mean calculation
    sum: i64,
    /// Running sum of squares for variance
    sum_sq: i64,
    /// Minimum observed
    min: i32,
    /// Maximum observed
    max: i32,
    /// Total samples seen (may exceed max_samples)
    total_seen: usize,
}

impl InsertSizeEstimator {
    pub fn new() -> Self {
        Self::with_capacity(10000)
    }

    pub fn with_capacity(max_samples: usize) -> Self {
        Self {
            samples: VecDeque::with_capacity(max_samples),
            max_samples,
            sum: 0,
            sum_sq: 0,
            min: i32::MAX,
            max: i32::MIN,
            total_seen: 0,
        }
    }

    /// Add an observed insert size
    pub fn add_sample(&mut self, insert_size: i32) {
        // Only track positive insert sizes
        if insert_size <= 0 {
            return;
        }

        self.total_seen += 1;
        self.min = self.min.min(insert_size);
        self.max = self.max.max(insert_size);

        let insert_i64 = insert_size as i64;
        self.sum += insert_i64;
        self.sum_sq += insert_i64 * insert_i64;

        // Keep sample buffer limited
        if self.samples.len() >= self.max_samples {
            // Remove oldest sample
            if let Some(old) = self.samples.pop_front() {
                let old_i64 = old as i64;
                self.sum -= old_i64;
                self.sum_sq -= old_i64 * old_i64;
            }
        }
        self.samples.push_back(insert_size);
    }

    /// Get current statistics
    pub fn stats(&self) -> InsertSizeStats {
        if self.samples.is_empty() {
            return InsertSizeStats::default();
        }

        let n = self.samples.len() as f64;
        let mean = self.sum as f64 / n;
        let variance = (self.sum_sq as f64 / n) - (mean * mean);
        let std_dev = variance.max(0.0).sqrt();

        // Calculate median
        let mut sorted: Vec<i32> = self.samples.iter().copied().collect();
        sorted.sort();
        let median = sorted[sorted.len() / 2];

        InsertSizeStats {
            mean,
            std_dev,
            min: self.min,
            max: self.max,
            count: self.total_seen,
            median,
        }
    }

    /// Check if we have enough samples for reliable estimation
    pub fn is_reliable(&self) -> bool {
        self.samples.len() >= 1000
    }

    /// Sample count
    pub fn sample_count(&self) -> usize {
        self.samples.len()
    }
}

impl Default for InsertSizeEstimator {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert_size_stats_default() {
        let stats = InsertSizeStats::default();
        assert_eq!(stats.mean, 300.0);
        assert_eq!(stats.std_dev, 50.0);
    }

    #[test]
    fn test_insert_size_expected() {
        let stats = InsertSizeStats {
            mean: 300.0,
            std_dev: 50.0,
            min: 100,
            max: 500,
            count: 1000,
            median: 300,
        };

        // Within 3 std dev: 150-450
        assert!(stats.is_expected(300));
        assert!(stats.is_expected(200));
        assert!(stats.is_expected(400));
        assert!(!stats.is_expected(100));
        assert!(!stats.is_expected(500));
    }

    #[test]
    fn test_insert_size_concordant() {
        let stats = InsertSizeStats {
            mean: 300.0,
            std_dev: 50.0,
            min: 100,
            max: 500,
            count: 1000,
            median: 300,
        };

        // Within 6 std dev: 0-600
        assert!(stats.is_concordant(300));
        assert!(stats.is_concordant(100));
        assert!(stats.is_concordant(500));
        assert!(!stats.is_concordant(700));
    }

    #[test]
    fn test_estimator_basic() {
        let mut estimator = InsertSizeEstimator::new();
        
        for _ in 0..100 {
            estimator.add_sample(300);
        }

        let stats = estimator.stats();
        assert!((stats.mean - 300.0).abs() < 0.01);
        assert!(stats.std_dev < 0.01); // All same value
    }

    #[test]
    fn test_estimator_distribution() {
        let mut estimator = InsertSizeEstimator::new();
        
        // Add samples centered around 300 with some spread
        for i in 200..400 {
            estimator.add_sample(i);
        }

        let stats = estimator.stats();
        assert!((stats.mean - 299.5).abs() < 1.0);
        assert!(stats.median >= 295 && stats.median <= 305);
        assert_eq!(stats.min, 200);
        assert_eq!(stats.max, 399);
    }

    #[test]
    fn test_estimator_ignores_negative() {
        let mut estimator = InsertSizeEstimator::new();
        
        estimator.add_sample(-100);
        estimator.add_sample(0);
        estimator.add_sample(300);

        assert_eq!(estimator.sample_count(), 1);
    }
}
