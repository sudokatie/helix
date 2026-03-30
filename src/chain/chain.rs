/// Seed chaining - connect collinear seeds into alignment chains

use super::Seed;

/// Configuration for seed chaining
#[derive(Debug, Clone)]
pub struct ChainConfig {
    /// Maximum gap in query between chained seeds
    pub max_query_gap: u32,
    /// Maximum gap in target between chained seeds
    pub max_target_gap: u32,
    /// Penalty per gap base
    pub gap_penalty: i32,
    /// Minimum chain score
    pub min_score: i32,
    /// Maximum diagonal difference for seeds to chain
    pub max_diagonal_diff: i64,
}

impl Default for ChainConfig {
    fn default() -> Self {
        Self {
            max_query_gap: 5000,
            max_target_gap: 5000,
            gap_penalty: 1,
            min_score: 20,
            max_diagonal_diff: 500,
        }
    }
}

/// A chain of collinear seeds
#[derive(Debug, Clone)]
pub struct Chain {
    /// Seeds in this chain (ordered by position)
    pub seeds: Vec<Seed>,
    /// Total chain score
    pub score: i32,
    /// Start position in query
    pub query_start: u32,
    /// End position in query
    pub query_end: u32,
    /// Start position in target
    pub target_start: u32,
    /// End position in target
    pub target_end: u32,
}

impl Chain {
    /// Create a new chain from a single seed
    #[allow(dead_code)]
    fn from_seed(seed: Seed) -> Self {
        Self {
            query_start: seed.query_pos,
            query_end: seed.query_pos + seed.length,
            target_start: seed.target_pos,
            target_end: seed.target_pos + seed.length,
            score: seed.length as i32,
            seeds: vec![seed],
        }
    }

    /// Extend chain with a new seed
    #[allow(dead_code)]
    fn extend(&mut self, seed: Seed, gap_cost: i32) {
        self.query_end = seed.query_pos + seed.length;
        self.target_end = seed.target_pos + seed.length;
        self.score += seed.length as i32 - gap_cost;
        self.seeds.push(seed);
    }

    /// Coverage of the query
    pub fn query_coverage(&self) -> u32 {
        self.query_end - self.query_start
    }

    /// Coverage of the target
    pub fn target_coverage(&self) -> u32 {
        self.target_end - self.target_start
    }
}

/// Chain seeds using dynamic programming
///
/// Seeds must be sorted by (target_pos, query_pos)
pub fn chain_seeds(seeds: &[Seed], config: &ChainConfig) -> Vec<Chain> {
    if seeds.is_empty() {
        return Vec::new();
    }

    let n = seeds.len();

    // DP: f[i] = best chain score ending at seed i
    let mut f: Vec<i32> = vec![0; n];
    // Backtrack: prev[i] = previous seed in best chain ending at i
    let mut prev: Vec<Option<usize>> = vec![None; n];

    // Initialize
    for i in 0..n {
        f[i] = seeds[i].length as i32;
    }

    // DP
    for i in 1..n {
        let si = &seeds[i];

        // Look back at previous seeds
        // In practice, we'd use a range tree for O(log n) lookups
        // For now, use simple O(n^2) with early termination
        let max_lookback = 1000.min(i);

        for j in (i.saturating_sub(max_lookback)..i).rev() {
            let sj = &seeds[j];

            // Check if seeds can be chained
            if let Some(score) = chain_score(sj, si, config) {
                let candidate = f[j] + score;
                if candidate > f[i] {
                    f[i] = candidate;
                    prev[i] = Some(j);
                }
            }
        }
    }

    // Backtrack to find chains
    extract_chains(&seeds, &f, &prev, config)
}

/// Calculate score for chaining seed j to seed i
fn chain_score(sj: &Seed, si: &Seed, config: &ChainConfig) -> Option<i32> {
    // Seeds must be in order (i comes after j)
    let query_gap = si.query_pos.checked_sub(sj.query_pos + sj.length)?;
    let target_gap = si.target_pos.checked_sub(sj.target_pos + sj.length)?;

    // Check maximum gaps
    if query_gap > config.max_query_gap || target_gap > config.max_target_gap {
        return None;
    }

    // Check diagonal difference (collinearity)
    let diag_diff = (si.diagonal() - sj.diagonal()).abs();
    if diag_diff > config.max_diagonal_diff {
        return None;
    }

    // Gap penalty based on the difference between query and target gaps
    // (indel-like gaps are penalized)
    let gap_diff = (query_gap as i64 - target_gap as i64).unsigned_abs() as u32;
    let gap_cost = (gap_diff as i32).saturating_mul(config.gap_penalty);

    // Score: seed length minus gap cost
    Some((si.length as i32).saturating_sub(gap_cost))
}

/// Extract chains from DP results using backtracking
fn extract_chains(
    seeds: &[Seed],
    f: &[i32],
    prev: &[Option<usize>],
    config: &ChainConfig,
) -> Vec<Chain> {
    let n = seeds.len();
    let mut used = vec![false; n];
    let mut chains = Vec::new();

    // Find chain endpoints (local maxima in f)
    let mut endpoints: Vec<(i32, usize)> = f.iter().copied().enumerate().map(|(i, s)| (s, i)).collect();
    endpoints.sort_by(|a, b| b.0.cmp(&a.0)); // Sort by score descending

    for (score, end_idx) in endpoints {
        if used[end_idx] || score < config.min_score {
            continue;
        }

        // Backtrack to find chain
        let mut chain_seeds = Vec::new();
        let mut idx = end_idx;

        loop {
            if used[idx] {
                break;
            }
            chain_seeds.push(seeds[idx]);
            used[idx] = true;

            match prev[idx] {
                Some(p) => idx = p,
                None => break,
            }
        }

        if chain_seeds.is_empty() {
            continue;
        }

        // Reverse to get correct order
        chain_seeds.reverse();

        // Build chain
        let first = chain_seeds[0];
        let last = chain_seeds.last().unwrap();

        let chain = Chain {
            query_start: first.query_pos,
            query_end: last.query_pos + last.length,
            target_start: first.target_pos,
            target_end: last.target_pos + last.length,
            score,
            seeds: chain_seeds,
        };

        chains.push(chain);
    }

    // Sort chains by score descending
    chains.sort_by(|a, b| b.score.cmp(&a.score));

    chains
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_chain_single_seed() {
        let seeds = vec![Seed::new(0, 100, 25)]; // Use 25 to exceed min_score of 20
        let config = ChainConfig::default();

        let chains = chain_seeds(&seeds, &config);

        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds.len(), 1);
        assert_eq!(chains[0].score, 25);
    }

    #[test]
    fn test_chain_collinear_seeds() {
        // Two seeds on same diagonal (collinear)
        let seeds = vec![
            Seed::new(0, 100, 15),
            Seed::new(20, 120, 15), // Same diagonal: 120-20 = 100-0 = 100
        ];
        let config = ChainConfig::default();

        let chains = chain_seeds(&seeds, &config);

        assert_eq!(chains.len(), 1);
        assert_eq!(chains[0].seeds.len(), 2);
        assert!(chains[0].score >= 30); // At least sum of seed lengths
    }

    #[test]
    fn test_chain_non_collinear() {
        // Two seeds on very different diagonals
        let seeds = vec![
            Seed::new(0, 100, 25),
            Seed::new(20, 1000, 25), // Different diagonal
        ];
        let config = ChainConfig {
            max_diagonal_diff: 100, // Restrictive
            min_score: 20,
            ..Default::default()
        };

        let chains = chain_seeds(&seeds, &config);

        // Should produce two separate chains
        assert_eq!(chains.len(), 2);
    }

    #[test]
    fn test_chain_gap_penalty() {
        // Seeds with a gap
        let seeds = vec![
            Seed::new(0, 100, 15),
            Seed::new(50, 155, 15), // 35 base gap in both, 5 base indel
        ];
        let config = ChainConfig {
            gap_penalty: 1,
            ..Default::default()
        };

        let chains = chain_seeds(&seeds, &config);

        assert_eq!(chains.len(), 1);
        // Score should be penalized for the gap difference
        assert!(chains[0].score < 30);
    }

    #[test]
    fn test_chain_min_score_filter() {
        let seeds = vec![Seed::new(0, 100, 10)];
        let config = ChainConfig {
            min_score: 20,
            ..Default::default()
        };

        let chains = chain_seeds(&seeds, &config);

        // Seed score (10) is below min_score (20)
        assert!(chains.is_empty());
    }

    #[test]
    fn test_chain_coverage() {
        let seeds = vec![
            Seed::new(0, 100, 15),
            Seed::new(20, 120, 15),
        ];
        let config = ChainConfig::default();

        let chains = chain_seeds(&seeds, &config);

        assert_eq!(chains[0].query_start, 0);
        assert_eq!(chains[0].query_end, 35); // 20 + 15
        assert_eq!(chains[0].query_coverage(), 35);
    }

    #[test]
    fn test_empty_seeds() {
        let seeds: Vec<Seed> = Vec::new();
        let config = ChainConfig::default();

        let chains = chain_seeds(&seeds, &config);

        assert!(chains.is_empty());
    }
}
