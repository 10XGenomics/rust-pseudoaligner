// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use debruijn::kmer;
use log::info;
use serde::Deserialize;

use anyhow::Error;
use bio::io::{fasta, fastq};
use docopt::Docopt;
use std::{env, fs};
use std::{path::PathBuf, str};

use debruijn_mapping::{
    build_index::build_index,
    mappability::{analyze_graph, write_mappability_tsv},
    pseudoaligner,
    pseudoaligner::process_reads,
    utils,
};

const PKG_NAME: &str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &str = env!("CARGO_PKG_VERSION");
const USAGE: &str = "
De-bruijn-mapping

Usage:
  pseudoaligner index [--num-threads=<n>] -i <index> <ref-fasta>
  pseudoaligner map [--num-threads=<n>] -i <index> <reads-fastq>
  pseudoaligner mappability [-o <outdir>] -i <index>
  pseudoaligner idxstats -i <index>
  pseudoaligner inspect -i <index> -c <counts> <genes>...
  pseudoaligner -h | --help | -v | --version

Options:
  -k --kmer-size      Kmer size to use - only 20 or 64 currently supported [default: 20].
  -n --num-threads N  Number of worker threads [default: 2]
  -o --outdir DIR     Output directory
  -h --help           Show this screen.
  -v --version        Show version.
";
// -l --long         Long output format (one line per read-transcript mapping)

#[derive(Clone, Debug, Deserialize)]
struct Args {
    arg_ref_fasta: String,
    arg_index: String,
    arg_reads_fastq: String,
    flag_outdir: Option<String>,
    flag_num_threads: usize,
    flag_kmer_size: usize,

    cmd_index: bool,

    cmd_map: bool,
    cmd_mappability: bool,
    cmd_idxstats: bool,

    // flag_long: bool,
    flag_version: bool,
    flag_v: bool,
}

enum KmerSetting {
    K20,
    K64,
}

fn main() -> Result<(), Error> {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());

    if args.flag_version || args.flag_v {
        println! {"{} {}", PKG_NAME, PKG_VERSION};
        return Ok(());
    }

    // initialize logger
    pretty_env_logger::init_timed();
    info!("Command line args:\n{:?}", args);

    let outdir = match &args.flag_outdir {
        Some(dir) => PathBuf::from(dir),
        None => env::current_dir()?,
    };
    fs::create_dir_all(&outdir)?;

    let km = match args.flag_kmer_size {
        20 => KmerSetting::K20,
        64 => KmerSetting::K64,
        v => {
            println!(
                "Kmer size = {} is not supported. Set kmer size to 20 or 64",
                v
            );
            return Ok(());
        }
    };

    if args.cmd_index {
        info!("Building index from fasta");
        let fasta = fasta::Reader::from_file(args.arg_ref_fasta)?;
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(fasta)?;

        match km {
            KmerSetting::K20 => {
                let index = build_index::<kmer::Kmer20>(
                    &seqs,
                    &tx_names,
                    &tx_gene_map,
                    args.flag_num_threads,
                )?;
                info!("Finished building index!");

                info!("Writing index to disk");
                utils::write_obj(&index, args.arg_index)?;
                info!("Finished writing index!");
            }
            KmerSetting::K64 => {
                let index = build_index::<kmer::Kmer64>(
                    &seqs,
                    &tx_names,
                    &tx_gene_map,
                    args.flag_num_threads,
                )?;
                info!("Finished building index!");

                info!("Writing index to disk");
                utils::write_obj(&index, args.arg_index)?;
                info!("Finished writing index!");
            }
        }
    } else if args.cmd_map {
        match km {
            KmerSetting::K20 => {
                info!("Reading index from disk");
                let index = utils::read_obj(args.arg_index)?;
                info!("Finished reading index!");

                info!("Mapping reads from fastq");
                let reads = fastq::Reader::from_file(args.arg_reads_fastq)?;
                process_reads::<kmer::Kmer20, _>(reads, &index, outdir, args.flag_num_threads)?;
            }
            KmerSetting::K64 => {
                info!("Reading index from disk");
                let index = utils::read_obj(args.arg_index)?;
                info!("Finished reading index!");

                info!("Mapping reads from fastq");
                let reads = fastq::Reader::from_file(args.arg_reads_fastq)?;
                process_reads::<kmer::Kmer64, _>(reads, &index, outdir, args.flag_num_threads)?;
            }
        }
    } else if args.cmd_mappability {
        match km {
            KmerSetting::K20 => {
                info!("Reading index from disk");
                let index = debruijn_mapping::utils::read_obj(args.arg_index)?;
                info!("Finished reading index!");
                info!("Analyzing de Bruijn graph");
                let records = analyze_graph::<kmer::Kmer20>(&index)?;
                info!("Finished analyzing!");
                info!("{} transcripts total", records.len());
                write_mappability_tsv(records, outdir)?;
            }
            KmerSetting::K64 => {
                info!("Reading index from disk");
                let index = debruijn_mapping::utils::read_obj(args.arg_index)?;
                info!("Finished reading index!");
                info!("Analyzing de Bruijn graph");
                let records = analyze_graph::<kmer::Kmer64>(&index)?;
                info!("Finished analyzing!");
                info!("{} transcripts total", records.len());
                write_mappability_tsv(records, outdir)?;
            }
        }
    } else if args.cmd_idxstats {
        match km {
            KmerSetting::K20 => {
                let index: pseudoaligner::Pseudoaligner<kmer::Kmer20> =
                    utils::read_obj(args.arg_index)?;

                use debruijn::Mer;

                for e in index.dbg.iter_nodes() {
                    let eqid = e.data();
                    let eq = &index.eq_classes[*eqid as usize];
                    println!("{}\t{}\t{}", e.node_id, e.sequence().len(), eq.len());
                }
            }
            KmerSetting::K64 => {
                let index: pseudoaligner::Pseudoaligner<kmer::Kmer64> =
                    utils::read_obj(args.arg_index)?;

                use debruijn::Mer;

                for e in index.dbg.iter_nodes() {
                    let eqid = e.data();
                    let eq = &index.eq_classes[*eqid as usize];
                    println!("{}\t{}\t{}", e.node_id, e.sequence().len(), eq.len());
                }
            }
        }
    }

    info!("Done!");
    Ok(())
}
