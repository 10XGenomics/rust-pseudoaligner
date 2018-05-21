extern crate debruijn;
extern crate bio;
extern crate clap;
extern crate itertools;
extern crate pdqsort;
extern crate boomphf;
extern crate smallvec;
extern crate pretty_env_logger;
extern crate bincode;
extern crate flate2;
extern crate failure;

#[macro_use]
extern crate log;

#[macro_use]
extern crate serde;

mod utils;

// Import some modules
use std::io;
use std::str;
use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::convert::AsRef;

use clap::{Arg, App};
use smallvec::SmallVec;
use bio::io::{fasta, fastq};

use debruijn::dna_string::*;
use debruijn::filter::filter_kmers;
use debruijn::graph::{DebruijnGraph};
use debruijn::{Dir, Kmer, Exts, kmer, Vmer};
use debruijn::compression::compress_kmers_with_hash;

const MIN_KMERS: usize = 1;
const STRANDED: bool = true;
const MEM_SIZE: usize = 1;
const REPORT_ALL_KMER: bool = false;
pub type PrimDataType = u32;
pub type KmerType = kmer::Kmer32;
pub type DataType = SmallVec<[PrimDataType; 4]>;


fn read_fasta(reader: fasta::Reader<File>)
              -> (DebruijnGraph<KmerType, DataType>,
                  boomphf::BoomHashMap2<KmerType, Exts, DataType>) {

    let summarizer = debruijn::filter::CountFilterSmallInt::new(MIN_KMERS);
    //let summarizer = filter::CountFilterSet::new(MIN_KMERS);
    let mut seqs = Vec::new();
    let mut trancript_counter: PrimDataType = 0;

    info!("Starting Reading the Fasta file\n");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        let dna_string = DnaString::from_dna_string( str::from_utf8(record.seq()).unwrap() );

        // obtain sequence and push into the relevant vector
        seqs.push((dna_string, Exts::empty(), trancript_counter));

        trancript_counter += 1;
        if trancript_counter % 10000 == 0 {
            eprint!("\r Done Reading {} sequences", trancript_counter);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        // looking for two transcripts
        // println!("{:?}", record.id());
        // if trancript_counter == 2 { break; }
    }
    eprintln!("");

    info!("Starting kmer filtering");
    //let (valid_kmers, obs_kmers): (Vec<(KmerType, (Exts, _))>, _) =
    //    filter::filter_kmers::<KmerType, _, _, _, _>(&seqs, summarizer, STRANDED);
    let (index, _) : (boomphf::BoomHashMap2<KmerType, Exts, _>, _) =
        filter_kmers::<KmerType, _, _, _, _>(&seqs, summarizer, STRANDED,
                                             REPORT_ALL_KMER, MEM_SIZE);

    //println!("Kmers observed: {}, kmers accepted: {}", obs_kmers.len(), valid_kmers.len());
    info!("Starting uncompressed de-bruijn graph construction");

    //println!("{:?}", index);

    let dbg = compress_kmers_with_hash(STRANDED, debruijn::compression::ScmapCompress::new(), &index).finish();
    info!("Done de-bruijn graph construction; ");

    let is_cmp = dbg.is_compressed();
    if is_cmp.is_some() {
        warn!("not compressed: nodes: {:?}", is_cmp);
        //dbg.print();
    }

    info!("Finished Indexing !");

    (dbg, index)
    //// TODO Should be added to ![cfg(test)] but doing here right now
    ////dbg.print_with_data();
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

fn process_reads(index: boomphf::BoomHashMap2<KmerType, Exts, DataType>,
                 dbg: DebruijnGraph<KmerType, DataType>,
                 reader: fastq::Reader<File>){

    let mut reads_counter = 0;
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        reads_counter += 1;

        let seqs = DnaString::from_dna_string( str::from_utf8(record.seq()).unwrap() );

        let mut eq_class: Vec<PrimDataType> = Vec::new();
        for kmer in seqs.iter_kmers() {
            //let (nid, _, _) = match dbg.find_link(kmer, Dir::Right){
            //    Some(links) => links,
            //    None => (std::usize::MAX, Dir::Right, false),
            //};
            //if nid != std::usize::MAX {
            //    let labels = dbg.get_node(nid).data();
            //    eq_class.extend(labels.clone().iter());
            //    pdqsort::sort(&mut eq_class);
            //    eq_class.dedup();
            //}
            let maybe_data = index.get(&kmer);
            match maybe_data {
                Some((_, ref labels)) => {
                    eq_class.extend(labels.clone().iter());
                    pdqsort::sort(&mut eq_class);
                    eq_class.dedup();
                },
                None => (),
            }
        }

        if reads_counter % 100000 == 0 {
            eprint!("\rDone Mapping {} reads", reads_counter);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        //println!("{:?} -> {:?}", record.id(), eq_class);
    }
    eprintln!();
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
             .help("Txome/Genome Input Fasta file")
             .required(true))
        .arg(Arg::with_name("reads")
             .short("r")
             .long("reads")
             .value_name("FILE")
             .help("Input Read Fastq file")
             .required(true))
        .arg(Arg::with_name("index")
             .short("i")
             .long("index")
             .value_name("FILE")
             .help("Index of the reference")
             .required(true))
        .arg(Arg::with_name("make")
             .help("tells to make the index")
             .short("m")
             .long("make")
             .requires("index"))
        .get_matches();
    pretty_env_logger::init();

    // Gets a value for config if supplied by user
    let fasta_file = matches.value_of("fasta").unwrap();
    info!("Path for reference FASTA: {}", fasta_file);

    let mut dbg: DebruijnGraph<KmerType, DataType>;
    let mut index: boomphf::BoomHashMap2<KmerType, Exts, DataType>;

    // obtain reader or fail with error (via the unwrap method)
    let index_file = matches.values_of("index").unwrap().collect();
    if matches.is_present("make") {
        let input_dump = utils::read_obj(index_file);
        match input_dump {
            Ok((dbg, index)) => (),
            Err(error) => (),
        }
    }
    else{
        let reader = fasta::Reader::from_file(fasta_file).unwrap();
        let output_dump = &index;
        utils::write_obj(output_dump, index_file);
    }

    // obtain reader or fail with error (via the unwrap method)
    let reads_file = matches.value_of("reads").unwrap();
    info!("Path for Reads FASTQ: {}\n\n", reads_file);

    let reads = fastq::Reader::from_file(reads_file).unwrap();
    process_reads(index, dbg, reads);
}
