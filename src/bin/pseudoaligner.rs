// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bio;
extern crate crossbeam_utils;
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
use debruijn_mapping::{bam, em};
                    


const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
De-bruijn-mapping

Usage:
  pseudoaligner index -i <index> <ref-fasta>
  pseudoaligner hla-index -i <index> <hla-fasta>
  pseudoaligner comb-index -i <index> <ref-fasta> <hla-fasta>
  pseudoaligner map -i <index> <reads-fastq>
  pseudoaligner map-bam [-o <outdir>]  -i <index> <bam> [<locus>]
  pseudoaligner hla-map     [-o <outdir>] -i <index> <reads-fastq>
  pseudoaligner mappability [-o <outdir>] -i <index>
  pseudoaligner idxstats -i <index>
  pseudoaligner em -i <index> -c <counts>
  pseudoaligner inspect -i <index> -c <counts> <genes>...
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
    arg_hla_fasta: String,
    arg_index: String,
    arg_counts: String,
    arg_bam: String,
    arg_locus: Option<String>,
    arg_genes: Vec<String>,
    arg_reads_fastq: String,
    flag_outdir: Option<String>,

    cmd_index: bool,
    cmd_hla_index: bool,
    cmd_comb_index: bool,

    cmd_em: bool,
    cmd_map: bool,
    cmd_map_bam: bool,
    cmd_mappability: bool,
    cmd_idxstats: bool,
    cmd_inspect: bool,

    // flag_long: bool,
    flag_version: bool,
    flag_v: bool,
}

const HLA_PATTERN: &'static str = "^HLA-.*";

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
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(fasta, None)?;
        let index = build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map
        )?;
        info!("Finished building index!");

        info!("Writing index to disk");
        utils::write_obj(&index, args.arg_index)?;
        info!("Finished writing index!");
    } else if args.cmd_hla_index {

        info!("Building index from fasta");
        let fasta = fasta::Reader::from_file(args.arg_hla_fasta)?;
        let (seqs, tx_names, tx_allele_map) = hla::read_hla_cds(fasta)?;

        let tx_gene_map = HashMap::new();

        let index = build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map
        )?;
        info!("Finished building index!");

        info!("Writing index to disk");
        utils::write_obj(&(index, tx_allele_map), args.arg_index)?;
        info!("Finished writing index!");
    } else if args.cmd_comb_index {

        info!("Building index from fasta");
        let fasta = fasta::Reader::from_file(args.arg_ref_fasta)?;

        // Load normal txs, with HLA genes filtered out.        
        let (mut seqs, mut tx_names, tx_gene_map) = utils::read_transcripts(fasta, Some(HLA_PATTERN))?;

        info!("Building index from fasta");
        let fasta = fasta::Reader::from_file(args.arg_hla_fasta)?;
        let (hla_seqs, hla_tx_names, hla_tx_allele_map) = hla::read_hla_cds(fasta)?;
        let tx_gene_map = HashMap::new();

        // WIP!!

        // Tack on HLA to 
        seqs.extend(hla_seqs);
        tx_names.extend(hla_tx_names);

        let index = build_index::<config::KmerType>(
            &seqs, &tx_names, &tx_gene_map
        )?;

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
    } else if args.cmd_em {


        let index: pseudoaligner::Pseudoaligner<config::KmerType> = utils::read_obj(args.arg_index)?;
        let mut eq_counts: bam::EqClassDb = utils::read_obj(&args.arg_counts)?;

        let eqclass_counts = eq_counts.eq_class_counts(index.tx_names.len());

        let weights = em::squarem(&eqclass_counts);

        use std::fs::File;
        use std::io::BufWriter;
        use std::io::Write;
        let mut weights_file = BufWriter::new(File::create("weights.tsv")?);

        let mut weight_names: Vec<(f64, &String, u32)> = 
            weights.
            into_iter().
            enumerate().
            map(|(i,w)| (w, &index.tx_names[i], i as u32)).
            collect();

        weight_names.sort_by(|(wa,_, _), (wb, _, _)| (-wa).partial_cmp(&-wb).unwrap());

        for (w, name, txid) in weight_names {
            if w > 3e-3 {
                // Find the tx's that overlap most strongly with this tx
                let mut tx_counts: HashMap<u32, u32> = HashMap::new();

                for (class, count) in &eqclass_counts.counts {
                    if class.binary_search(&txid).is_ok() {
                        for tx in class {
                            let c = tx_counts.entry(*tx).or_insert(0);
                            *c += count;
                        }
                    }
                }

                let my_counts = *(tx_counts.iter().filter(|(ttx, _)| *ttx == &txid).next().unwrap().1);

                let mut neighbor_counts: Vec<_> = tx_counts.into_iter().filter(|(ttx, _)| *ttx != txid).collect();
                neighbor_counts.sort_by_key(|(_,c)| -(*c as i32));


                println!("X: {} {}", name, my_counts);
                for i in 0 .. std::cmp::min(12, neighbor_counts.len()) {
                    let (neighbor_id, count) = neighbor_counts[i].clone();
                    if (count as f64) < 0.5 * my_counts as f64 { break }
                    println!("_: {}  {}", &index.tx_names[neighbor_id as usize], count);
                }
            }

            writeln!(weights_file, "{}\t{}", name, w)?;
        }

    } else if args.cmd_inspect {

        let index: pseudoaligner::Pseudoaligner<config::KmerType> = utils::read_obj(args.arg_index)?;
        let mut eq_counts: bam::EqClassDb = utils::read_obj(&args.arg_counts)?;

        let eqclass_counts = eq_counts.eq_class_counts(index.tx_names.len());

        let mut name_map = HashMap::new();
        for (i, name) in index.tx_names.iter().enumerate() {
            name_map.insert(name, i);
        }

        println!("genes: {:?}", args.arg_genes);
        for s in args.arg_genes {
            let id = name_map.get(&s).expect("couldn't find gene").clone() as u32;

            println!("===== {} =====", s);
            for (class, count) in &eqclass_counts.counts {
                if count > &10 {
                    if class.binary_search(&id).is_ok() {
                        for tx in class {
                            println!("g: {}", index.tx_names[*tx as usize]);
                        }
                        println!("count: {}", count);
                    }
                }
            }
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