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
use std::collections::HashMap;
use bio::io::{fasta, fastq};
use docopt::Docopt;
use failure::Error;
use std::{env, fs};
use std::{path::PathBuf, str};

use debruijn_mapping::{config, utils, hla};
use debruijn_mapping::{build_index::build_index,
                       pseudoaligner,
                       pseudoaligner::process_reads,
                       mappability::analyze_graph};
use debruijn_mapping::bam;
                    


const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
De-bruijn-mapping

Usage:
  pseudoaligner index -i <index> <ref-fasta>
  pseudoaligner hla-index -i <index> <ref-fasta>
  pseudoaligner map -i <index> <reads-fastq>
  pseudoaligner map-bam [-o <outdir>]  -i <index> <locus> <bam>
  pseudoaligner hla-map     [-o <outdir>] -i <index> <reads-fastq>
  pseudoaligner mappability [-o <outdir>] -i <index>
  pseudoaligner idxstats -i <index>
  pseudoaligner -h | --help | -v | --version

Options:
  -o --outdir DIR   Output directory
  -h --help         Show this screen.
  -v --version         Show version.
";
  // -l --long         Long output format (one line per read-transcript mapping)


#[derive(Clone, Debug, Deserialize)]
struct Args {
    arg_ref_fasta: String,
    arg_index: String,
    arg_bam: String,
    arg_locus: String,
    arg_reads_fastq: String,
    flag_outdir: Option<String>,
    cmd_index: bool,
    cmd_hla_index: bool,
    cmd_map: bool,
    cmd_map_bam: bool,
    cmd_mappability: bool,
    cmd_idxstats: bool,

    // flag_long: bool,
    flag_version: bool,
    flag_v: bool,
}

fn main() -> Result<(), Error> {
    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.deserialize())
                            .unwrap_or_else(|e| e.exit());

    if args.flag_version || args.flag_v {
        println!{"{} {}", PKG_NAME, PKG_VERSION};
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

    if args.cmd_index {
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
    } else if args.cmd_hla_index {

        info!("Building index from fasta");
        let fasta = fasta::Reader::from_file(args.arg_ref_fasta)?;
        let (seqs, tx_names, tx_allele_map) = hla::read_hla_cds(fasta)?;

        let tx_gene_map = HashMap::new();

        let index = build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map
        )?;
        info!("Finished building index!");

        info!("Writing index to disk");
        utils::write_obj(&(index, tx_allele_map), args.arg_index)?;
        info!("Finished writing index!");

    } else if args.cmd_map {
        info!("Reading index from disk");
        let index = utils::read_obj(args.arg_index)?;
        info!("Finished reading index!");

        info!("Mapping reads from fastq");
        let reads = fastq::Reader::from_file(args.arg_reads_fastq)?;
        process_reads::<config::KmerType, _>(reads, &index, outdir)?;
    } else if args.cmd_map_bam {
        map_bam(&args.clone())?;
    } else if args.cmd_mappability {
        info!("Reading index from disk");
        let index = debruijn_mapping::utils::read_obj(args.arg_index)?;
        info!("Finished reading index!");
        info!("Analyzing de Bruijn graph");
        let records = analyze_graph::<config::KmerType>(&index)?;
        info!("Finished analyzing!");
        info!("{} transcripts total", records.len());
        utils::write_mappability_tsv(records, outdir)?;
    } else if args.cmd_idxstats {

        let index: pseudoaligner::Pseudoaligner<config::KmerType> = utils::read_obj(args.arg_index)?;

        use debruijn::Mer;

        for e in index.dbg.iter_nodes() {
            let eqid = e.data();
            let eq = &index.eq_classes[*eqid as usize];
            println!("{}\t{}\t{}", e.node_id, e.sequence().len(), eq.len());
        }
    }

    info!("Done!");
    Ok(())
}


fn map_bam(args: &Args) -> Result<(), Error> {

    info!("Reading index from disk");
    let index = utils::read_obj(&args.arg_index)?;
    info!("Finished reading index!");

    info!("Mapping reads from BAM");

    let path = PathBuf::from(args.flag_outdir.as_ref().unwrap());
    //let p2: PathBuf = args.flag_outdir.as_ref().unwrap().into();
    bam::map_bam(&args.arg_bam, index, &args.arg_locus, &path)?;

    Ok(())
}