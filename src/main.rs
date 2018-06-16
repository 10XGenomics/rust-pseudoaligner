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
extern crate crossbeam;

#[macro_use]
extern crate log;

#[macro_use]
extern crate serde;

mod utils;
mod docks;
mod work_queue;
mod config;

// Import some modules
use std::io;
use std::str;
use std::fs::File;
use std::ops::Add;
use std::io::Write;
use std::fmt::Debug;
use std::hash::Hash;

use std::sync::mpsc;
use std::sync::{Arc, Mutex};

use num::One;
use num::Num;
use num::NumCast;
use clap::{Arg, App};
use serde::Serialize;
use bio::io::{fasta, fastq};
use serde::de::DeserializeOwned;

use debruijn::{Exts};
use debruijn::dna_string::*;

use config::{U8_MAX, U16_MAX, U32_MAX, MAX_WORKER, BUCKET_SIZE_THRESHOLD, MIN_KMERS};
use config::{DocksUhs, KmerType};

fn read_fasta(reader: fasta::Reader<File>)
              -> Vec<Vec<DnaString>> {
    let mut seqs = Vec::new();
    let mut transcript_counter = 0;

    info!("Starting Reading the Fasta file\n");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        let dna_string = DnaString::from_dna_only_string( str::from_utf8(record.seq()).unwrap() );

        // obtain sequence and push into the relevant vector
        seqs.push(dna_string);

        transcript_counter += 1;
        if transcript_counter % 100 == 0 {
            print!("\r Done Reading {} sequences", transcript_counter);
            io::stdout().flush().ok().expect("Could not flush stdout");
        }
        // looking for two transcripts
        // println!("{:?}", record.id());
        // if trancript_counter == 2 { break; }
        //warn!("Reading only one Chromosome");
        //if transcript_counter > 1 { break };
    }
    println!();
    info!("Done Reading the Fasta file; Found {} sequences", transcript_counter);

    seqs
}

fn filter_kmers_callback(seqs: Vec<Vec<DnaString>>, index_file: &str,
                         uhs: DocksUhs) {
    let seqs_len = seqs.len();

    // Based on the number of sequences chose the right primary datatype
    match seqs_len {
        1 ...U8_MAX  => {
            info!("Using 8 bit variable for storing the data.");
            call_filter_kmers(&seqs, index_file,
                              &uhs, u8::min_value());
        },
        U8_MAX ... U16_MAX => {
            info!("Using 16 bit variable for storing the data.");
            call_filter_kmers(&seqs, index_file,
                              &uhs, u16::min_value());
        },
        U16_MAX ... U32_MAX => {
            info!("Using 32 bit variable for storing the data.");
            call_filter_kmers::<u32>(&seqs, index_file,
                                     &uhs, u32::min_value());
        },
        _ => {
            error!("Too many ({}) sequneces to handle.", seqs_len);
        },
    };
}

fn call_filter_kmers<S>(seqs: &Vec<Vec<DnaString>>, index_file: &str,
                        uhs: &DocksUhs, _seq_id: S)
where S: Clone + Hash + Eq + Debug + Ord + Serialize + One + Add<Output=S>
    + Send + Sync + Num + NumCast + DeserializeOwned {
    info!("Starting Bucketing");
    let (tx, rx) = mpsc::channel();
    let num_seqs = seqs.len();
    let queue = Arc::new(work_queue::WorkQueue::new());
    let mut buckets: Vec<Vec<(DnaStringSlice, Exts, S)>> = vec![Vec::new(); uhs.len()];

    info!("Spawning {} threads for Bucketing.", MAX_WORKER);
    crossbeam::scope(|scope| {

        for _ in 0 .. MAX_WORKER {
            let tx = tx.clone();
            let queue = Arc::clone(&queue);

            scope.spawn(move || {
                loop {
                    let done = queue.len();
                    if done % 10 == 0 {
                        print!("\rDone Bucketing {}% of the reference sequences",
                               std::cmp::min(100, done*100/num_seqs));
                        io::stdout().flush().ok().expect("Could not flush stdout");
                    }

                    // If work is available, do that work.
                    match queue.get_work(seqs) {
                        Some((seq, head)) => {
                            let thread_data =  work_queue::run(seq, uhs);
                            tx.send((thread_data, head)).expect("Could not send data!");
                          },
                        None => { break; },
                    };
                } // end loop
            });
        }

        let mut missed_bases_counter: usize = 0;

        for _ in 0..num_seqs {
            let (seq_data, head) = rx.recv().unwrap();
            let (bucket_slices, missed_bases) = seq_data;

            let bit_head: S = num::cast(head).unwrap();
            missed_bases_counter += missed_bases;
            for (bucket_id, slices) in bucket_slices {
                buckets[bucket_id as usize].push((slices, Exts::empty(), bit_head.clone()));
            }
        }

        println!();
        info!("Bucketing successfully finished.");
        warn!("Missed total {} bases", missed_bases_counter);
    });

    // separating small and big buckets since call to Boomphf new is expensive
    let mut big_buckets: Vec<Vec<(DnaStringSlice, Exts, S)>> = Vec::new();
    let mut small_bucket: Vec<(DnaStringSlice, Exts, S)> = Vec::new();
    for bucket in buckets {
        let num_elem = bucket.len();
        if num_elem != 0 {
            if num_elem > BUCKET_SIZE_THRESHOLD {
                big_buckets.push(bucket);
            }
            else{
                for elem in bucket {
                    small_bucket.push(elem);
                }
            }
        }
    }

    // do all small buckets at once
    let mut num_buckets = big_buckets.len();
    big_buckets.insert(num_buckets/2, small_bucket);
    num_buckets += 1;

    let mut dbgs = Vec::new();
    let summarizer = Arc::new(debruijn::filter::CountFilterEqClass::new(MIN_KMERS));
    {
        info!("Starting per-bucket De-bruijn Graph Creation");
        let (tx, rx) = mpsc::channel();

        let atomic_buckets = Arc::new(Mutex::new(big_buckets));
        let queue = Arc::new(work_queue::WorkQueue::new());

        info!("Spawning {} threads for Analyzing.", MAX_WORKER);
        crossbeam::scope(|scope| {
            for _ in 0 .. MAX_WORKER {
                let tx = tx.clone();
                let queue = Arc::clone(&queue);
                let bucket_ref = Arc::clone(&atomic_buckets);
                let summarizer_ref = Arc::clone(&summarizer);

                scope.spawn(move || {
                    loop {
                        let done = queue.len();
                        if done % 10 == 0 {
                            print!("\rDone Analyzing {}% of the buckets",
                                   done*100/num_buckets);
                            io::stdout().flush().ok().expect("Could not flush stdout");
                        }

                        // If work is available, do that work.
                        match queue.get_rev_work(&bucket_ref) {
                            Some((bucket_data, _)) => {
                                let thread_data = work_queue::analyze(bucket_data, &summarizer_ref);
                                tx.send(thread_data).expect("Could not send data!");
                            },
                            None => { break; },
                        };
                    } // end loop
                });
            }
        });

        for _ in 0..num_buckets {
            let dbg = rx.recv().unwrap();
            dbgs.push(dbg.unwrap());
        }
    }

    println!();
    info!("Done seprate de-bruijn graph construction; ");
    info!("Starting merge");
    //println!("{:?}", summarizer);

    let full_dbg = work_queue::merge_graphs(dbgs);
    let eq_classes = Arc::try_unwrap(summarizer).unwrap().get_eq_classes();

    utils::Index::dump(full_dbg,
                       eq_classes.into_inner().expect("Can't unwrap eqclass"),
                       index_file);
}

fn process_reads<S>(index: utils::Index<KmerType, S>,
                    reader: fastq::Reader<File>)
where S: Clone + Ord + PartialEq + Debug + Sync + Send + Hash + Serialize + DeserializeOwned {
    info!("Done Reading index");
    info!("Starting Multi-threaded Mapping");

    let (tx, rx) = mpsc::channel();
    let atomic_reader = Arc::new(Mutex::new(reader.records()));

    let dbg = index.get_dbg();
    let left_phf = index.get_left_phf();
    let right_phf = index.get_right_phf();
    let eq_classes = index.get_eq_classes();

    info!("Spawning {} threads for Mapping.", MAX_WORKER);
    crossbeam::scope(|scope| {
        for _ in 0 .. MAX_WORKER {
            let tx = tx.clone();
            let reader = Arc::clone(&atomic_reader);

            scope.spawn(move || {
                loop {
                    // If work is available, do that work.
                    match work_queue::get_next_record(&reader) {
                        Some(result_record) => {
                            let record = match result_record {
                                Ok(record) => record,
                                Err(err) => panic!("Error {:?} in reading fastq", err),
                            };

                            let seqs = DnaString::from_dna_string( str::from_utf8(record.seq()).unwrap() );
                            let eq_class = work_queue::map(seqs, dbg, eq_classes);

                            tx.send(Some(eq_class)).expect("Could not send data!");
                        },
                        None => { break; },
                    }; //end-match
                } // end loop

                // send None to tell receiver that the thread exited
                tx.send(None).expect("Could not send data!");
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
                },
                Some(eq_class) => {
                    read_counter += 1;
                    if eq_class.len() > 0 {
                        mapped_read_counter += 1;
                    }

                    if read_counter % 100000 == 0 {
                        let frac_mapped = mapped_read_counter as f32 * 100.0 / read_counter as f32;
                        print!("\rDone Mapping {} reads w/ Rate: {}",
                               read_counter, frac_mapped);
                        io::stdout().flush().ok().expect("Could not flush stdout");
                    }

                    eprintln!("{:?}", eq_class);
                } // end-Some
            } // end-match
        } // end-for
    }); //end crossbeam

    println!();
    info!("Done Mapping Reads");
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

    // obtain reader or fail with error (via the unwrap method)
    let index_file = matches.values_of("index").unwrap().next().unwrap();

    if matches.is_present("make") {
        warn!("Creating the index, can take little time.");
        let uhs = docks::read_uhs();

        // if index not found then create a new one
        let reader = fasta::Reader::from_file(fasta_file).unwrap();
        let seqs = read_fasta(reader);

        //Set up the filter_kmer call based on the number of sequences.
        filter_kmers_callback(seqs, index_file, uhs);

        info!("Finished Indexing !");
    }
    else{
        // import the index if already present.
        info!("Reading index from File: {:?}", index_file);
        let data_type_file_name = index_file.to_owned() + "/type.bin";
        let data_type: usize = utils::read_obj(data_type_file_name).expect("Can't read data type");

        // obtain reader or fail with error (via the unwrap method)
        let reads_file = matches.value_of("reads").unwrap();
        info!("Path for Reads FASTQ: {}\n\n", reads_file);
        let reads = fastq::Reader::from_file(reads_file).unwrap();

        match data_type {
            1 => {
                info!("Read u8 for eqclass data type");
                let index = utils::Index::<KmerType, u8>::read(index_file);
                process_reads(index, reads);
            },
            2 => {
                info!("Read u16 for eqclass data type");
                let index = utils::Index::<KmerType, u16>::read(index_file);
                process_reads(index, reads);
            },
            4 => {
                info!("Read u32 for eqclass data type");
                let index = utils::Index::<KmerType, u32>::read(index_file);
                process_reads(index, reads);
            },
            _ => panic!("read unidentified data type with size => {:?}", data_type),
        };
    }
    info!("Finished Processing !")
}

#[cfg(test)]
mod tests{
    use std;
    use utils;
    use bincode;
    use smallvec::SmallVec;
    use debruijn::{Dir, Kmer, Exts, kmer};

    pub type KmerType = kmer::Kmer32;

    #[test]
    fn test_kmer_search() {
        let index_file = "/mnt/home/avi.srivastava/rust_avi/rust-utils-10x/sc_mapping/unit_test/test.small.index";
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
