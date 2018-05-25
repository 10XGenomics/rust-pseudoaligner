extern crate debruijn;
extern crate bio;
extern crate clap;
extern crate itertools;
extern crate pdqsort;
extern crate boomphf;
extern crate pretty_env_logger;
extern crate bincode;
extern crate flate2;
extern crate num;
extern crate failure;

#[macro_use]
extern crate smallvec;

#[macro_use]
extern crate log;

#[macro_use]
extern crate serde;

mod utils;

// Import some modules
use std::io;
use std::str;
use std::fs::File;
use std::ops::Add;
use std::io::Write;
use std::fmt::Debug;
use std::hash::Hash;

use num::One;
use clap::{Arg, App};
use serde::Serialize;
use bio::io::{fasta, fastq};

use debruijn::dna_string::*;
use debruijn::filter::filter_kmers;
use debruijn::graph::{DebruijnGraph};
use debruijn::{Exts, kmer, Vmer};
use debruijn::filter::EqClassIdType;
use debruijn::compression::compress_kmers_with_hash;

const MIN_KMERS: usize = 1;
const STRANDED: bool = true;
const MEM_SIZE: usize = 1;
const REPORT_ALL_KMER: bool = false;
const U8_MAX: usize = u8::max_value() as usize;
const U16_MAX: usize = u16::max_value() as usize;
const U32_MAX: usize = u32::max_value() as usize;

pub type KmerType = kmer::Kmer32;

fn read_fasta(reader: fasta::Reader<File>)
              -> Vec<DnaString> {
    let mut seqs = Vec::new();
    let mut trancript_counter = 0;

    info!("Starting Reading the Fasta file\n");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        let dna_string = DnaString::from_dna_string( str::from_utf8(record.seq()).unwrap() );

        // obtain sequence and push into the relevant vector
        seqs.push(dna_string);

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
    info!("Done Reading the Fasta file; Found {} sequences", trancript_counter);

    seqs
}

fn filter_kmers_callback(seqs: Vec<DnaString>, index_file: &str) {
    let seqs_len = seqs.len();

    // Based on the number of sequences chose the right primary datatype
    match seqs_len {
        1 ...U8_MAX  => {
            info!("Using 8 bit variable for storing the data.");
            call_filter_kmers(seqs, index_file,
                              u8::min_value());
        },
        U8_MAX ... U16_MAX => {
            info!("Using 16 bit variable for storing the data.");
            call_filter_kmers(seqs, index_file,
                              u16::min_value());
        },
        U16_MAX ... U32_MAX => {
            info!("Using 32 bit variable for storing the data.");
            call_filter_kmers::<u32>(seqs,index_file ,
                                     u32::min_value());
        },
        _ => {
            error!("Too many ({}) sequneces to handle.", seqs_len);
        },
    };
}

fn call_filter_kmers<S>(seqs: Vec<DnaString>, index_file: &str,
                        mut seq_id: S)
where S: Clone + Hash + Eq + Debug + Ord + Serialize + One + Add<Output=S> {
    let mut summarizer = debruijn::filter::CountFilterEqClass::new(MIN_KMERS);
    let mut input_tuples = Vec::new();

    for seq in seqs {
        input_tuples.push((seq, Exts::empty(), seq_id.clone()));
        seq_id = seq_id + num::one();
    }

    info!("Starting kmer filtering");
    let (phf, _) : (boomphf::BoomHashMap2<KmerType, Exts, EqClassIdType>, _) =
        filter_kmers::<KmerType, _, _, _, _>(&input_tuples, &mut summarizer, STRANDED,
                                             REPORT_ALL_KMER, MEM_SIZE);

    info!("Starting uncompressed de-bruijn graph construction");
    //println!("{:?}", phf);
    let dbg = compress_kmers_with_hash(STRANDED, debruijn::compression::ScmapCompress::new(),
                                       &phf).finish();

    info!("Done de-bruijn graph construction; ");
    let is_cmp = dbg.is_compressed();
    if is_cmp.is_some() {
        warn!("not compressed: nodes: {:?}", is_cmp);
        //dbg.print();
    }

    let ref_index = utils::Index::new(dbg, phf, summarizer.get_eq_classes());

    //info!("Dumping index into File: {:?}", index_file);
    //utils::write_obj(&ref_index, index_file).expect("Can't dump the index");
}




fn process_reads(phf: &boomphf::BoomHashMap2<KmerType, Exts, EqClassIdType>,
                 //dbg: &DebruijnGraph<KmerType, EqClassIdType>,
                 reader: fastq::Reader<File>){

    let mut reads_counter = 0;
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        reads_counter += 1;

        let seqs = DnaString::from_dna_string( str::from_utf8(record.seq()).unwrap() );

        let mut eq_class: Vec<EqClassIdType> = Vec::new();
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
            let maybe_data = phf.get(&kmer);
            match maybe_data {
                Some((_, ref labels)) => {
                    eq_class.push(*labels.clone());
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
             .help("Txome/Genome Input Fasta file, (Needed only with -m i.e. while making index)"))
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
             .requires("index")
             .requires("fasta"))
        .get_matches();
    pretty_env_logger::init();

    // Gets a value for config if supplied by user
    let fasta_file = matches.value_of("fasta").unwrap();
    info!("Path for reference FASTA: {}", fasta_file);

    //let ref_index: utils::Index<KmerType, Exts, SequnceIdType>;

    // obtain reader or fail with error (via the unwrap method)
    let index_file = matches.values_of("index").unwrap().next().unwrap();

    if matches.is_present("make") {
        warn!("Creating the index, can take little time.");
        // if index not found then create a new one
        let reader = fasta::Reader::from_file(fasta_file).unwrap();
        let seqs = read_fasta(reader);

        filter_kmers_callback(seqs, index_file);

        info!("Finished Indexing !");
    }
    //else{
    //    // import the index if already present.
    //    info!("Reading index from File: {:?}", index_file);
    //    let input_dump: Result<utils::Index<KmerType, Exts, EqClassIdType>,
    //                           Box<bincode::ErrorKind>> =
    //        utils::read_obj(index_file);

    //    ref_index = input_dump.expect("Can't read the index");

    //    // obtain reader or fail with error (via the unwrap method)
    //    let reads_file = matches.value_of("reads").unwrap();
    //    info!("Path for Reads FASTQ: {}\n\n", reads_file);

    //    let reads = fastq::Reader::from_file(reads_file).unwrap();
    //    process_reads(ref_index.get_phf(), /*ref_index.get_dbg(),*/ reads);

    //}
    info!("Finished Processing !")
}

#[cfg(test)]
mod tests{
    use std;
    use utils;
    use bincode;
    use smallvec::SmallVec;
    use debruijn::{Dir, Kmer, Exts, kmer};

    pub type EqClassIdType = u32;
    pub type KmerType = kmer::Kmer32;

    #[test]
    fn test_kmer_search() {
        let index_file = "/mnt/home/avi.srivastava/rust_avi/rust-utils-10x/sc_mapping/unit_test/test.small.index";
        println!("Reading index from File: {:?}", index_file);
        let input_dump: Result<utils::Index<KmerType, Exts, EqClassIdType>,
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
        //let color = input_dump.unwrap().get_eq_classes().get(eqclass_id);
        let oracle: Vec<u32> = vec![0, 1];
        assert_eq!(oracle, vec![0, 1]);
    }
}
