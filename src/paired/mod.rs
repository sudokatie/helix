// Paired-end read support
// Handles paired alignment, insert size estimation, and concordant pair detection

mod pair;
mod insert;

pub use pair::{ReadPair, PairedAligner, PairOrientation};
pub use insert::{InsertSizeEstimator, InsertSizeStats};
