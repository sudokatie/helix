use clap::{Parser, Subcommand};
use helix::chain::{chain_seeds, extend_chain, find_seeds, ChainConfig, ExtendConfig, SeedConfig};
use helix::index::{save_index, IndexBuilder, IndexConfig, load_index};
use helix::output::{SamHeader, SamRecord, SamWriter, FLAG_REVERSE};
use helix::seq::{open_fasta, open_fastq, decode_sequence, reverse_complement};
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
            let loaded_index = load_index(&index_path)?;
            
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

            // Add reference sequences from index
            for ref_info in &loaded_index.references {
                // Extract reference name (handle empty names)
                let name = if ref_info.name.is_empty() {
                    "unnamed".to_string()
                } else {
                    ref_info.name.clone()
                };
                header.add_reference(name, ref_info.length);
            }

            sam_writer.write_header(&header)?;

            // Process reads
            let reader = open_fastq(&reads)?;
            let mut read_count = 0;
            let mut aligned_count = 0;

            // Alignment configuration
            let seed_config = SeedConfig::default();
            let chain_config = ChainConfig::default();
            let extend_config = ExtendConfig::default();
            let min_score_threshold = 30;

            for result in reader {
                let seq = result?;
                let seq_decoded = decode_sequence(&seq.seq);
                let seq_str = String::from_utf8_lossy(&seq_decoded).to_string();
                let qual_str = String::from_utf8_lossy(&seq.qual).to_string();

                // Try forward orientation
                let mut best_alignment = try_align(
                    &seq.seq,
                    &loaded_index,
                    &seed_config,
                    &chain_config,
                    &extend_config,
                );
                let mut is_reverse = false;

                // Try reverse complement if forward alignment is poor
                let rc_seq = reverse_complement(&seq.seq);
                let rc_alignment = try_align(
                    &rc_seq,
                    &loaded_index,
                    &seed_config,
                    &chain_config,
                    &extend_config,
                );

                // Pick the better orientation
                if let Some(ref rc_aln) = rc_alignment {
                    if best_alignment.as_ref().map_or(true, |fwd| rc_aln.score > fwd.score) {
                        best_alignment = rc_alignment;
                        is_reverse = true;
                    }
                }

                // Output record
                let record = if let Some(aln) = best_alignment {
                    if aln.score >= min_score_threshold {
                        aligned_count += 1;

                        // Compute MAPQ
                        let mapq = aln.mapq;

                        let mut flag: u16 = 0;
                        if is_reverse {
                            flag |= FLAG_REVERSE;
                        }

                        SamRecord {
                            qname: seq.id.clone(),
                            flag,
                            rname: aln.ref_name.clone(),
                            pos: (aln.ref_start + 1) as u32, // 1-based SAM position
                            mapq,
                            cigar: aln.cigar.to_string(),
                            rnext: "*".to_string(),
                            pnext: 0,
                            tlen: 0,
                            seq: seq_str,
                            qual: qual_str,
                            tags: vec![
                                ("AS".to_string(), format!("i:{}", aln.score)),
                            ],
                        }
                    } else {
                        SamRecord::new_unmapped(seq.id.clone(), seq_str, qual_str)
                    }
                } else {
                    SamRecord::new_unmapped(seq.id.clone(), seq_str, qual_str)
                };
                sam_writer.write_record(&record)?;

                read_count += 1;
                if read_count % 10000 == 0 {
                    eprintln!("  processed {} reads...", read_count);
                }
            }

            eprintln!("Aligned {} / {} reads ({:.1}%)", aligned_count, read_count,
                      100.0 * aligned_count as f64 / read_count.max(1) as f64);
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

/// Try to align a sequence to the index
fn try_align(
    seq: &[u8],
    index: &helix::index::Index,
    seed_config: &SeedConfig,
    chain_config: &ChainConfig,
    extend_config: &ExtendConfig,
) -> Option<helix::chain::Alignment> {
    // Step 1: Find seeds
    let seeds = find_seeds(seq, index, seed_config);
    if seeds.is_empty() {
        return None;
    }

    // Step 2: Chain seeds
    let chains = chain_seeds(&seeds, chain_config);
    if chains.is_empty() {
        return None;
    }

    // Step 3: Extend best chain
    let best_chain = &chains[0]; // chains are sorted by score descending
    let mut alignment = extend_chain(best_chain, seq, index, extend_config)?;

    // Step 4: Compute MAPQ from best and second-best scores
    let best_score = alignment.score;
    let second_best = chains.get(1).map(|c| c.score).unwrap_or(0);
    alignment.mapq = compute_mapq(best_score, second_best);

    Some(alignment)
}

/// Compute MAPQ from alignment score and suboptimal score
fn compute_mapq(best_score: i32, subopt_score: i32) -> u8 {
    if best_score <= 0 {
        return 0;
    }

    // MAPQ = min(60, -10 * log10(1 - best_score / (best_score + subopt_score)))
    // Simplified: higher difference between best and second-best = higher MAPQ
    let total = best_score + subopt_score;
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
