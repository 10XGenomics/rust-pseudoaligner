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
extern crate csv;
extern crate rayon;

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
use std::collections::HashMap;

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
use config::{READ_COVERAGE_THRESHOLD};
use config::{DocksUhs, KmerType};

fn read_fasta(reader: fasta::Reader<File>,
              tgmap_file: &str)
              -> (Vec<Vec<DnaString>>, Vec<usize>, Vec<String>) {
    let mut gene_id: usize = 0;
    let tgmap_reader = File::open(tgmap_file).expect("can't read tgmap file");
    let mut tgmap: HashMap<String, usize> = HashMap::new();
    let mut gene_id_map: HashMap<String, usize> = HashMap::new();

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(tgmap_reader);

    for maybe_line in rdr.records() {
        let curr_gene_id: usize;
        let records = maybe_line.expect("tgmap read error");

        let transcript_name = &records[0];
        let gene_name = &records[1];

        if gene_id_map.contains_key( gene_name ) {
            curr_gene_id = gene_id_map[gene_name];
        }
        else{
            curr_gene_id = gene_id.clone();
            gene_id_map.insert(gene_name.to_string(), gene_id);
            gene_id += 1;
        }

        tgmap.insert(transcript_name.to_string(), curr_gene_id);
    }

    let mut seqs = Vec::new();
    let mut transcript_counter = 0;
    let mut tgmap_vec = Vec::new();

    info!("Starting Reading the Fasta file\n");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result.unwrap();
        let dna_string = DnaString::from_dna_only_string( str::from_utf8(record.seq()).unwrap() );
        let record_id: Vec<&str> = record.id().split('|').collect();

        let gene_id = tgmap.get(record_id[0]).expect("can't find fasta entry in tgMap");
        tgmap_vec.push(gene_id.clone());

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

    let mut gene_id_vec = vec!["".to_string(); gene_id_map.len()];
    for (gene, gid) in gene_id_map {
        gene_id_vec[ gid ] = gene;
    }

    (seqs, tgmap_vec, gene_id_vec)
}

fn filter_kmers_callback(seqs: Vec<Vec<DnaString>>, index_file: &str,
                         uhs: DocksUhs, tgmap: Vec<usize>,
                         gene_order: Vec<String>) {
    let seqs_len = seqs.len();

    // Based on the number of sequences chose the right primary datatype
    match seqs_len {
        1 ...U8_MAX  => {
            info!("Using 8 bit variable for storing the data.");
            call_filter_kmers(&seqs, index_file,
                              &uhs, u8::min_value(),
                              tgmap, gene_order);
        },
        U8_MAX ... U16_MAX => {
            info!("Using 16 bit variable for storing the data.");
            call_filter_kmers(&seqs, index_file,
                              &uhs, u16::min_value(),
                              tgmap, gene_order);
        },
        U16_MAX ... U32_MAX => {
            info!("Using 32 bit variable for storing the data.");
            call_filter_kmers::<u32>(&seqs, index_file,
                                     &uhs, u32::min_value(),
                                     tgmap, gene_order);
        },
        _ => {
            error!("Too many ({}) sequneces to handle.", seqs_len);
        },
    };
}

fn call_filter_kmers<S>(seqs: &Vec<Vec<DnaString>>, index_file: &str,
                        uhs: &DocksUhs, _seq_id: S, tgmap: Vec<usize>,
                        gene_order: Vec<String>)
where S: Clone + Hash + Eq + Debug + Ord + Serialize + One + Add<Output=S>
    + Send + Sync + Num + NumCast + DeserializeOwned {
    info!("Starting Bucketing");
    let (tx, rx) = mpsc::sync_channel(MAX_WORKER);
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
                    // If work is available, do that work.
                    match queue.get_work(seqs) {
                        Some((seq, head)) => {
                            let thread_data =  work_queue::run(seq, uhs);
                            tx.send((thread_data, head)).expect("Could not send data!");

                            let done = queue.len();
                            if done % 10 == 0 {
                                print!("\rDone Bucketing {}% of the reference sequences",
                                       std::cmp::min(100, done*100/num_seqs));
                                io::stdout().flush().ok().expect("Could not flush stdout");
                            }
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
            let gene_id = tgmap.get(head).expect("transcript id out of range");

            let bit_head: S = num::cast(gene_id.clone()).unwrap();
            missed_bases_counter += missed_bases;
            for (bucket_id, slices, exts) in bucket_slices {
                buckets[bucket_id as usize].push((slices, exts, bit_head.clone()));
            }
        }

        println!();
        info!("Bucketing successfully finished.");
        warn!("Missed total {} bases", missed_bases_counter);
    }); //end-crossbeam

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
        let (tx, rx) = mpsc::sync_channel(MAX_WORKER);

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
                        // If work is available, do that work.
                        match queue.get_rev_work(&bucket_ref) {
                            Some((bucket_data, _)) => {
                                let thread_data = work_queue::analyze(bucket_data, &summarizer_ref);
                                tx.send(thread_data).expect("Could not send data!");

                                let done = queue.len();
                                if done % 10 == 0 {
                                    print!("\rDone Analyzing {}% of the buckets",
                                           done*100/num_buckets);
                                    io::stdout().flush().ok().expect("Could not flush stdout");
                                }
                            },
                            None => { break; },
                        };
                    } // end loop
                });
            }

            for _ in 0..num_buckets {
                let dbg = rx.recv().unwrap();
                dbgs.push(dbg.unwrap());
            }
        }); // end-crossbeam
    }

    println!();
    info!("Done seprate de-bruijn graph construction; ");
    info!("Starting merge");

    // Thread pool Configuration for calling BOOMphf
    rayon::ThreadPoolBuilder::new().num_threads(MAX_WORKER).build_global().unwrap();

    //println!("{:?}", summarizer);
    let full_dbg = work_queue::merge_graphs(dbgs);
    let eq_classes = Arc::try_unwrap(summarizer).ok().unwrap().get_eq_classes();

    utils::Index::dump(full_dbg, gene_order,
                       eq_classes,
                       index_file);
}

fn process_reads<S>(index: utils::Index<KmerType, S>,
                    reader: fastq::Reader<File>)
where S: Clone + Ord + PartialEq + Debug + Sync + Send + Hash + Serialize + DeserializeOwned {
    info!("Done Reading index");
    info!("Starting Multi-threaded Mapping");

    let (tx, rx) = mpsc::sync_channel(MAX_WORKER);
    let atomic_reader = Arc::new(Mutex::new(reader.records()));

    let dbg = index.get_dbg();
    let phf = index.get_phf();
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

                            let dna_string = str::from_utf8(record.seq()).unwrap();
                            let seqs = DnaString::from_dna_string( dna_string );
                            let read_data = work_queue::map(seqs, dbg, eq_classes, phf);

                            let wrapped_read_data = match read_data {
                                Some((eq_class, coverage)) => {
                                    Some((record.id().to_owned(), eq_class, coverage))
                                },
                                None => Some(("".to_string(), Vec::new(), 0)),
                            };

                            tx.send(wrapped_read_data).expect("Could not send data!");
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
                Some(read_data) => {
                    read_counter += 1;
                    if read_counter % 100000 == 0 {
                        let frac_mapped = mapped_read_counter as f32 * 100.0 / read_counter as f32;
                        print!("\rDone Mapping {} reads w/ Rate: {}",
                               read_counter, frac_mapped);
                        io::stdout().flush().ok().expect("Could not flush stdout");
                    }

                    let (read_id, eq_class, coverage) = read_data;
                    if coverage >= READ_COVERAGE_THRESHOLD && eq_class.len() > 0 {
                        mapped_read_counter += 1;
                        eprintln!("{}=>{:?},{}", read_id, eq_class, coverage);
                    }
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
             .help("Input Read Fastq file"))
        .arg(Arg::with_name("tgMap")
             .short("t")
             .long("tgMap")
             .value_name("FILE")
             .help("Transcript to Gene Mapping file")
             .requires("fasta"))
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
             .requires("tgMap")
             .requires("index")
             .requires("fasta"))
        .get_matches();
    pretty_env_logger::init();

    // obtain reader or fail with error (via the unwrap method)
    let index_file = matches.values_of("index").unwrap().next().unwrap();

    if matches.is_present("make") {
        warn!("Creating the index, can take little time.");
        let uhs = docks::read_uhs();

        // Gets a value for config if supplied by user
        let fasta_file = matches.value_of("fasta").unwrap();
        info!("Path for reference FASTA: {}", fasta_file);

        // obtain reader or fail with error (via the unwrap method)
        let tgmap_file = matches.values_of("tgMap").unwrap().next().expect("no tgmap file found");

        // if index not found then create a new one
        let reader = fasta::Reader::from_file(fasta_file).unwrap();
        let (seqs, tgmap, gene_order) = read_fasta(reader, tgmap_file);

        //Set up the filter_kmer call based on the number of sequences.
        filter_kmers_callback(seqs, index_file, uhs, tgmap, gene_order);

        info!("Finished Indexing !");
    }
    else{
        // import the index if already present.
        info!("Reading index from File: {:?}", index_file);
        let data_type = utils::get_data_type(index_file);

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
