use clap::{Parser, Subcommand};
use helix::index::{save_index, IndexBuilder, IndexConfig, load_index};
use helix::output::{SamHeader, SamRecord, SamWriter};
use helix::seq::{open_fasta, open_fastq, decode_sequence};
use helix::util::detect_simd;
use std::fs::File;
use std::io::{BufWriter, stdout};
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

        /// Minimizer window size
        #[arg(short, long, default_value = "10")]
        window: usize,
    },

    /// Align reads to a reference
    Align {
        /// Index prefix
        index: PathBuf,

        /// Reads FASTQ file
        reads: PathBuf,

        /// Output SAM file (- for stdout)
        #[arg(short, long, default_value = "-")]
        output: String,

        /// Number of threads
        #[arg(short, long, default_value = "1")]
        threads: usize,
    },

    /// Show version and detected SIMD level
    Info,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index {
            reference,
            output,
            kmer,
            window,
        } => {
            eprintln!("Building index from {:?}", reference);
            eprintln!("  k-mer size: {}", kmer);
            eprintln!("  window size: {}", window);

            let config = IndexConfig::new(kmer, window);
            let mut builder = IndexBuilder::new(config);

            // Read references
            let reader = open_fasta(&reference)?;
            let mut ref_count = 0;
            let mut total_bases = 0u64;

            for result in reader {
                let reference = result?;
                total_bases += reference.seq.len() as u64;
                builder.add_reference(&reference);
                ref_count += 1;
                
                if ref_count % 100 == 0 {
                    eprintln!("  processed {} references...", ref_count);
                }
            }

            eprintln!("  {} references, {} bases", ref_count, total_bases);

            let index = builder.build();
            eprintln!("  {} unique minimizers", index.hash.len());

            // Save index
            let index_path = output.with_extension("hlix");
            save_index(&index, &index_path)?;
            eprintln!("Index saved to {:?}", index_path);
        }

        Commands::Align {
            index,
            reads,
            output,
            threads: _threads,
        } => {
            eprintln!("Loading index from {:?}", index);
            let index_path = index.with_extension("hlix");
            let _index = load_index(&index_path)?;
            
            eprintln!("Aligning reads from {:?}", reads);

            // Set up output
            let writer: Box<dyn std::io::Write> = if output == "-" {
                Box::new(stdout())
            } else {
                Box::new(BufWriter::new(File::create(&output)?))
            };
            let mut sam_writer = SamWriter::new(writer);

            // Write header
            let mut header = SamHeader::new();
            header.program = Some("helix".to_string());
            // TODO: Add reference sequences from index
            sam_writer.write_header(&header)?;

            // Process reads
            let reader = open_fastq(&reads)?;
            let mut read_count = 0;

            for result in reader {
                let seq = result?;
                
                // TODO: Actually align the read
                // For now, output as unmapped
                let record = SamRecord::new_unmapped(
                    seq.id.clone(),
                    String::from_utf8_lossy(&decode_sequence(&seq.seq)).to_string(),
                    String::from_utf8_lossy(&seq.qual).to_string(),
                );
                sam_writer.write_record(&record)?;
                
                read_count += 1;
                if read_count % 10000 == 0 {
                    eprintln!("  processed {} reads...", read_count);
                }
            }

            eprintln!("Aligned {} reads", read_count);
        }

        Commands::Info => {
            println!("helix {}", env!("CARGO_PKG_VERSION"));
            println!("SIMD level: {}", detect_simd().name());
            
            #[cfg(target_arch = "x86_64")]
            {
                if is_x86_feature_detected!("avx2") {
                    println!("  AVX2: available");
                } else {
                    println!("  AVX2: not available");
                }
                if is_x86_feature_detected!("avx512f") {
                    println!("  AVX-512: available");
                } else {
                    println!("  AVX-512: not available");
                }
            }
            
            #[cfg(target_arch = "aarch64")]
            {
                println!("  NEON: available (always on aarch64)");
            }
        }
    }

    Ok(())
}
