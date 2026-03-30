pub mod seq;
pub mod index;
pub mod align;
pub mod chain;
pub mod output;
pub mod util;

use thiserror::Error;

#[derive(Debug, Error)]
pub enum HelixError {
    #[error("Invalid FASTQ format: {0}")]
    FastqFormat(String),

    #[error("Invalid FASTA format: {0}")]
    FastaFormat(String),

    #[error("Index not found: {0}")]
    IndexNotFound(std::path::PathBuf),

    #[error("Unsupported SIMD: this binary requires AVX2")]
    UnsupportedSimd,

    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
}

pub type Result<T> = std::result::Result<T, HelixError>;
