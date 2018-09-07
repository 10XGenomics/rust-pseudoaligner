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
use std::str;

use bio::io::{fasta, fastq};

use failure::Error;

use docopt::Docopt;

use debruijn_mapping::{build_index, config, utils};
use debruijn_mapping::pseudoaligner::{process_reads};

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
De-bruijn-mapping

Usage:
  pseudoaligner index <output> <ref-fasta>
  pseudoaligner map <index> <reads-fastq>
  pseudoaligner -h | --help
  pseudoaligner --version 

Options:
  -h --help     Show this screen.
  --version     Show version.
";

#[derive(Debug, Deserialize)]
struct Args {
    arg_output: String,
    arg_ref_fasta: String,
    arg_index: String,
    arg_reads_fastq: String,
    cmd_index: bool,
    cmd_map: bool,
    flag_version: bool,
}

fn main() -> Result<(), Error> {
    let args: Args = Docopt::new(USAGE)
                            .and_then(|d| d.deserialize())
                            .unwrap_or_else(|e| e.exit());

    // initializing logger
    pretty_env_logger::init_timed();

    info!("Command line args:\n{:?}", args);

    if args.flag_version {
        println!{"{} {}", PKG_NAME, PKG_VERSION};
    } else if args.cmd_index {
        let fasta_fn = args.arg_ref_fasta;
        let output_index_fn = args.arg_output;

        warn!("Creating the index, can take little time.");
        info!("Path for reference FASTA: {}", fasta_fn);

        // if index not found then create a new one
        let reader = fasta::Reader::from_file(fasta_fn)?;
        let (seqs, _tx_ids, _tx_gene_map) = utils::read_fasta(reader)?;

        //Set up the filter_kmer call based on the number of sequences.
        let index = build_index::build_pseudoaligner_index::<config::KmerType>(&seqs);
        utils::write_obj(&index, output_index_fn)?;

        info!("Finished Indexing !");
    } else if args.cmd_map {
        // import the index
        let index_fn = args.arg_index;
        let input_reads_fn = args.arg_reads_fastq;
        
        info!("Reading index from File: {:?}", index_fn);
        let index = utils::read_obj(index_fn)?;

        info!("Path for Reads FASTQ: {}\n\n", input_reads_fn);
        let reads = fastq::Reader::from_file(input_reads_fn)?;

        process_reads::<config::KmerType>(&index, reads);
    }

    info!("Finished Processing !");
    Ok(())
}

/*
#[cfg(test)]
mod tests{
    use std;
    use utils;
    use bincode;
    use debruijn::{Dir, Kmer, Exts, kmer};

    pub type KmerType = kmer::Kmer32;

    const 

    #[test]
    fn test_kmer_search() {
        let tx_file = "/mnt/home/avi.srivastava/rust_avi/rust-utils-10x/sc_mapping/unit_test/test.small.index";
        println!("Reading index from File: {:?}", index_file);
        let input_dump: Result<utils::Index<KmerType, Exts, u8>,
                               Box<bincode::ErrorKind>> =
            utils::read_obj(index_file);

        let ref_index = input_dump.expect("Can't read the index");

        println!("Starting Unit test for color extraction");
        let test_kmer = KmerType::from_ascii(b"GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTT");
        let (nid, _, _) = match ref_index.get_dbg().find_link(test_kmer, Dir::Right){
            Some(links) => links,
            None => (std::usize::MAX, Dir::Right, false),
        };
        if nid == std::usize::MAX {
            eprintln!("ERROR");
        }
        println!("Found Colors are");
        let eqclass_id = ref_index.get_dbg().get_node(nid).data();
        let eq_classes = ref_index.get_eq_classes();
        let labels = eq_classes[*eqclass_id as usize].clone();
        assert_eq!(labels, vec![0, 1]);
    }
}
*/
