extern crate debruijn;
extern crate bio;
extern crate clap;
extern crate itertools;
extern crate pdqsort;
extern crate boomphf;
extern crate fxhash;
extern crate bit_set;
extern crate heapsize;
extern crate num_traits;
extern crate smallvec;

// helper functions for this project
mod utils;

// Import some modules
use std::io;
use std::str;
use std::fs::File;
use std::io::Write;

use clap::{Arg, App};
use bio::io::fasta;

use debruijn::{filter};
use debruijn::dna_string::*;
use debruijn::{Dir, Kmer, Exts, kmer};
use debruijn::compression::{compress_kmers};

pub type KmerType = kmer::Kmer32;
const MIN_KMERS: usize = 1;
const STRANDED: bool = true;
const REPORT_ALL_KMER: bool = false;

fn read_fasta(reader: fasta::Reader<File>) -> () {

    let summarizer = utils::CountFilterSmallInt::new(MIN_KMERS);
    //let summarizer = filter::CountFilterSet::new(MIN_KMERS);
    let mut seqs = Vec::new();
    let mut trancript_counter: u32 = 0;

    println!("Starting Reading the Fasta file");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        let dna_string = DnaString::from_dna_string( str::from_utf8(record.seq()).unwrap() );

        // obtain sequence and push into the relevant vector
        seqs.push((dna_string, Exts::empty(), trancript_counter));

        trancript_counter += 1;
        if trancript_counter % 10000 == 0 {
            print!("\r Done Reading {} transcripts", trancript_counter);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        // looking for two transcripts
        // println!("{:?}", record.id());
        // if trancript_counter == 2 { break; }
    }

    println!("\nStarting kmer filtering");
    let index: utils::BoomHashMap<KmerType, Exts, _> =
        utils::filter_kmers_with_mphf::<KmerType, _, _, _, _>(seqs, summarizer, STRANDED,
                                                              REPORT_ALL_KMER, 1);

    //println!("Kmers observed: {}, kmers accepted: {}", obs_kmers.len(), valid_kmers.len());
    println!("Starting uncompressed de-bruijn graph construction");

    //println!("{:?}", valid_kmers);

    let dbg = utils::compress_kmers_with_mphf(STRANDED, utils::ScmapCompress::new(), &index).finish();
    println!("Done de-bruijn graph construction; ");

    let is_cmp = dbg.is_compressed();
    if is_cmp.is_some() {
        println!("not compressed: nodes: {:?}", is_cmp);
        //dbg.print();
    }

    println!("Finished Indexing !");

    //// TODO Should be added to ![cfg(test)] but doing here right now
    //dbg.print_with_data();
    //println!("Starting Unit test for color extraction");
    //let test_kmer = KmerType::from_ascii(b"GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTT");
    //let (nid, _, _) = match dbg.find_link(test_kmer, Dir::Right){
    //    Some(links) => links,
    //    None => (std::usize::MAX, Dir::Right, false),
    //};
    //if nid == std::usize::MAX {
    //    eprintln!("ERROR");
    //}
    //println!("Found Colors are");
    //println!("{:?}", dbg.get_node(nid).data());
}

fn main() {
    let matches = App::new("De-bruijn-mapping")
        .version("1.0")
        .author("Avi S. <avi.srivastava@10xgenomics.com>")
        .about("De-bruijn graph based lightweight mapping for single-cell data")
        .arg(Arg::with_name("fasta")
             .short("f")
             .long("fasta")
             .value_name("FILE")
             .help("Genome Input file")
             .required(true))
        .get_matches();

    // Gets a value for config if supplied by user
    let fasta_file = matches.value_of("fasta").unwrap();
    println!("Path for FASTA: {}", fasta_file);

    // obtain reader or fail with error (via the unwrap method)
    let reader = fasta::Reader::from_file(fasta_file).unwrap();

    read_fasta(reader);
}
