use clap::{Parser, Subcommand};
use std::path::PathBuf;

#[derive(Parser)]
#[command(name = "helix")]
#[command(about = "DNA sequence aligner with SIMD acceleration")]
#[command(version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build an index from a reference FASTA file
    Index {
        /// Reference FASTA file
        reference: PathBuf,

        /// Output index prefix
        #[arg(short, long)]
        output: PathBuf,

        /// K-mer size for indexing
        #[arg(short, long, default_value = "15")]
        kmer: usize,
    },

    /// Align reads to a reference
    Align {
        /// Index prefix
        index: PathBuf,

        /// Reads FASTQ file
        reads: PathBuf,

        /// Output SAM file
        #[arg(short, long)]
        output: PathBuf,

        /// Number of threads
        #[arg(short, long, default_value = "1")]
        threads: usize,
    },
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index {
            reference,
            output,
            kmer,
        } => {
            println!(
                "Building index from {:?} with k={} -> {:?}",
                reference, kmer, output
            );
            // TODO: Implement indexing
        }
        Commands::Align {
            index,
            reads,
            output,
            threads,
        } => {
            println!(
                "Aligning {:?} to index {:?} with {} threads -> {:?}",
                reads, index, threads, output
            );
            // TODO: Implement alignment
        }
    }
}
