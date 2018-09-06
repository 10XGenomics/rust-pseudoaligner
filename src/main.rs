// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bincode;
extern crate bio;
extern crate boomphf;
extern crate crossbeam;
extern crate csv;
extern crate debruijn;
extern crate docopt;
extern crate failure;
extern crate flate2;
extern crate itertools;
extern crate num;
extern crate pretty_env_logger;
extern crate rayon;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate log;

#[macro_use]
extern crate serde;

mod build_index;
mod config;
mod pseudoaligner;
mod utils;

// Import some modules
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::io::Write;
use std::str;

use std::sync::mpsc;
use std::sync::{Arc, Mutex};

use bio::io::{fasta, fastq};

use debruijn::dna_string::*;
use debruijn::Kmer;

use pseudoaligner::Pseudoaligner;

use config::MAX_WORKER;
use config::READ_COVERAGE_THRESHOLD;

use failure::Error;
use flate2::read::MultiGzDecoder;
use std::io::{BufRead, BufReader};
use std::path::Path;

use docopt::Docopt;

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

/// Open a (possibly gzipped) file into a BufReader.
fn _open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<BufRead>, Error> {
    let r = File::open(p.as_ref())?;

    if p.as_ref().extension().unwrap() == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(32 * 1024, gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32 * 1024, r);
        Ok(Box::new(buf_reader))
    }
}

fn read_fasta(
    reader: fasta::Reader<File>,
) -> Result<(Vec<DnaString>, Vec<String>, HashMap<String, String>), Error> {
    let mut seqs = Vec::new();
    let mut transcript_counter = 0;
    let mut tx_ids = Vec::new();
    let mut tx_to_gene_map = HashMap::new();

    info!("Starting Reading the Fasta file\n");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result?;

        // Sequence
        let dna_string = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
        seqs.push(dna_string);

        let headers: Vec<&str> = record.id().split('|').collect();

        let tx_id = headers[0].to_string();
        let gene_id = headers[1].to_string();
        tx_ids.push(tx_id.clone());
        tx_to_gene_map.insert(tx_id, gene_id);

        transcript_counter += 1;
        if transcript_counter % 100 == 0 {
            print!("\r Done Reading {} sequences", transcript_counter);
            io::stdout().flush().expect("Could not flush stdout");
        }
    }

    println!();
    info!(
        "Done Reading the Fasta file; Found {} sequences",
        transcript_counter
    );

    Ok((seqs, tx_ids, tx_to_gene_map))
}

fn process_reads<K>(index: &Pseudoaligner<K>, reader: fastq::Reader<File>)
where
    K: Kmer + Sync + Send,
{
    info!("Done Reading index");
    info!("Starting Multi-threaded Mapping");

    let (tx, rx) = mpsc::sync_channel(MAX_WORKER);
    let atomic_reader = Arc::new(Mutex::new(reader.records()));

    info!("Spawning {} threads for Mapping.\n", MAX_WORKER);
    crossbeam::scope(|scope| {
        for _ in 0..MAX_WORKER {
            let tx = tx.clone();
            let reader = Arc::clone(&atomic_reader);

            scope.spawn(move || {
                loop {
                    // If work is available, do that work.
                    match get_next_record(&reader) {
                        Some(result_record) => {
                            let record = match result_record {
                                Ok(record) => record,
                                Err(err) => panic!("Error {:?} in reading fastq", err),
                            };

                            let dna_string = str::from_utf8(record.seq()).unwrap();
                            let seq = DnaString::from_dna_string(dna_string);
                            let read_data = index.map_read(&seq);

                            let wrapped_read_data = match read_data {
                                Some((eq_class, coverage)) => {
                                    if coverage >= READ_COVERAGE_THRESHOLD && eq_class.is_empty() {
                                        Some((true, record.id().to_owned(), eq_class, coverage))
                                    } else {
                                        Some((false, record.id().to_owned(), eq_class, coverage))
                                    }
                                }
                                None => Some((false, record.id().to_owned(), Vec::new(), 0)),
                            };

                            tx.send(wrapped_read_data).expect("Could not send data!");
                        }
                        None => {
                            // send None to tell receiver that the queue ended
                            tx.send(None).expect("Could not send data!");
                            break;
                        }
                    }; //end-match
                } // end loop
            }); //end-scope
        } // end-for

        let mut read_counter: usize = 0;
        let mut mapped_read_counter: usize = 0;
        let mut dead_thread_count = 0;

        for eq_class in rx.iter() {
            match eq_class {
                None => {
                    dead_thread_count += 1;
                    if dead_thread_count == MAX_WORKER {
                        drop(tx);
                        // can't continue with a flag check
                        // weird Rusty way !
                        // Consume whatever is remaining
                        // Not worrying about counters; hunch is their
                        // should be less
                        for eq_class in rx.iter() {
                            eq_class.map_or((), |eq_class| eprintln!("{:?}", eq_class));
                        }
                        break;
                    }
                }
                Some(read_data) => {
                    println!("{:?}", read_data);

                    if read_data.0 {
                        mapped_read_counter += 1;
                    }

                    read_counter += 1;
                    if read_counter % 1_000_000 == 0 {
                        let frac_mapped = mapped_read_counter as f32 * 100.0 / read_counter as f32;
                        eprint!(
                            "\rDone Mapping {} reads w/ Rate: {}",
                            read_counter, frac_mapped
                        );
                        io::stderr().flush().expect("Could not flush stdout");
                    }
                } // end-Some
            } // end-match
        } // end-for
    }); //end crossbeam

    eprintln!();
    info!("Done Mapping Reads");
}

pub fn get_next_record<R>(
    reader: &Arc<Mutex<fastq::Records<R>>>,
) -> Option<Result<fastq::Record, std::io::Error>>
where
    R: std::io::Read,
{
    let mut lock = reader.lock().unwrap();
    lock.next()
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
        let reader = fasta::Reader::from_file(fasta_fn).unwrap();
        let (seqs, _tx_ids, _tx_gene_map) = read_fasta(reader)?;

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
