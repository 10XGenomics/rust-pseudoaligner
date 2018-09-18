// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bio;
extern crate crossbeam;
extern crate debruijn;
extern crate debruijn_mapping;
extern crate docopt;
extern crate failure;
extern crate pretty_env_logger;
extern crate rayon;

#[macro_use]
extern crate log;

#[macro_use]
extern crate serde;

// Import some modules
use bio::io::{fasta, fastq};
use docopt::Docopt;
use failure::Error;
use std::str;

use debruijn_mapping::{config, utils};
use debruijn_mapping::{build_index::build_index, pseudoaligner::process_reads};

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
De-bruijn-mapping

Usage:
  pseudoaligner index -i <index> <ref-fasta>
  pseudoaligner map [options] -i <index> <reads-fastq>
  pseudoaligner -h | --help | --version

Options:
  -l --long         Long output format (one line per read-transcript mapping)
  -o --outdir DIR   Output directory
  -h --help         Show this screen.
  --version         Show version.
";

#[derive(Debug, Deserialize)]
struct Args {
    arg_ref_fasta: String,
    arg_index: String,
    arg_reads_fastq: String,
    flag_outdir: Option<String>,
    cmd_index: bool,
    cmd_map: bool,
    flag_long: bool,
    flag_version: bool,
}

fn main() -> Result<(), Error> {
    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.deserialize())
                            .unwrap_or_else(|e| e.exit());
    // initialize logger
    pretty_env_logger::init_timed();
    info!("Command line args:\n{:?}", args);

    if args.flag_version {
        println!{"{} {}", PKG_NAME, PKG_VERSION};
    } else if args.cmd_index {
        info!("Building index from fasta");
        let fasta = fasta::Reader::from_file(args.arg_ref_fasta)?;
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(fasta)?;
        let index = build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map
        )?;
        info!("Finished building index!");

        info!("Writing index to disk");
        utils::write_obj(&index, args.arg_index)?;
        info!("Finished writing index!");
    } else if args.cmd_map {
        info!("Reading index from disk");
        let index = utils::read_obj(args.arg_index)?;
        info!("Finished reading index!");

        info!("Mapping reads from fastq");
        let outdir = args.flag_outdir;
        let reads = fastq::Reader::from_file(args.arg_reads_fastq)?;
        process_reads::<config::KmerType, _>(reads, &index, outdir)?;
        info!("Finished mapping reads!");
    }

    info!("Done!");
    Ok(())
}
