// Part of code taken from
// https://gist.github.com/LeoTindall/e6d40782b05dc8ac40faf3a0405debd3
use std;
use std::hash::Hash;
use std::fmt::Debug;
use std::sync::{Mutex, Arc};
use std::sync::atomic::{AtomicUsize, Ordering};

use KmerType;
use STRANDED;
use REPORT_ALL_KMER;
use MEM_SIZE;

use boomphf;
use pdqsort;
use debruijn::*;
use debruijn::graph::*;
use debruijn::filter::*;
use bio::io::{fastq};
use debruijn::compression::*;
use docks::{L, DocksUhs, generate_msps};
use debruijn::dna_string::{DnaString, DnaStringSlice};

pub const MAX_WORKER: usize = 20;

pub struct WorkQueue{
    head : AtomicUsize,
}

impl WorkQueue{
    pub fn new() -> Self {
        Self {
            head: AtomicUsize::new(0),
        }
    }

    pub fn get_work<'a, T>(&self, seqs: &'a Vec<T>)
                           -> Option<(&'a T, usize)> {
        let old_head = self.head.fetch_add(1, Ordering::SeqCst);
        match seqs.get(old_head) {
            None => None,
            Some(contigs) => Some( (contigs, old_head) ),
        }
    }

    pub fn get_rev_work<T:Clone>(&self, atomic_seqs: &Arc<Mutex<Vec<T>>>)
                                  -> Option<(T, usize)> {
        let old_head = self.head.fetch_add(1, Ordering::SeqCst);
        let mut seqs = atomic_seqs.lock().expect("Can't unlock Buckets");
        match seqs.pop() {
            None => None,
            Some(contigs) => Some( (contigs, old_head) ),
        }
    }

    pub fn len(&self) -> usize {
        self.head.load(Ordering::SeqCst)
    }
}

pub fn run<'a>(contigs: &'a Vec<DnaString>, uhs: &DocksUhs)
               -> (std::vec::Vec<(u16, DnaStringSlice<'a>)>, usize){
    // One FASTA entry possibly broken into multiple contigs
    // based on the location of `N` int he sequence.
    let mut missed_bases_counter = 0;
    let mut bucket_slices = Vec::new();

    for seq in contigs {
        let seq_len = seq.len();
        if seq_len >= L {
            let msps = generate_msps( &seq, uhs );
            for msp in msps{
                let bucket_id = msp.bucket();
                if bucket_id > uhs.len() as u16{
                    panic!("Small bucket size");
                }

                let slice = seq.slice(msp.start(), msp.end());
                bucket_slices.push((bucket_id, slice));
            }
        }
        else{
            missed_bases_counter += seq_len;
        }
    }
    (bucket_slices, missed_bases_counter)
}

pub fn analyze<S:Clone+Eq+Hash+Ord+Debug>( bucket_data: Vec<(DnaStringSlice, Exts, S)>,
                                           summarizer: &Arc<CountFilterEqClass<S>>)
                                           -> Option<BaseGraph<KmerType, EqClassIdType>> {
    if bucket_data.len() > 0 {
        //println!("{:?}", bucket_data);
        // run filter_kmer
        let (phf, _) : (boomphf::BoomHashMap2<KmerType, Exts, EqClassIdType>, _) =
            filter_kmers::<KmerType, _, _, _, _>(&bucket_data, &summarizer, STRANDED,
                                                 REPORT_ALL_KMER, MEM_SIZE);

        //println!("{:?}", phf);
        // compress the graph
        let dbg = compress_kmers_with_hash(STRANDED, ScmapCompress::new(), phf);

        return Some(dbg);
    }
    None
}

pub fn merge_graphs( uncompressed_dbgs: Vec<BaseGraph<KmerType, EqClassIdType>> )
                     -> DebruijnGraph<KmerType, EqClassIdType> {
    //println!("{:?}", uncompressed_dbgs);
    // make a combined graph
    let combined_graph = BaseGraph::combine(uncompressed_dbgs.into_iter()).finish();

    //println!("{:?}", combined_graph);
    // compress the graph
    let dbg_graph = compress_graph(STRANDED, ScmapCompress::new(),
                                   combined_graph, None);

    //println!("{:?}", dbg_graph);
    dbg_graph
}

pub fn get_next_record<R>(reader: &Arc<Mutex<fastq::Records<R>>>)
                      -> Option<Result<fastq::Record, std::io::Error>>
where R: std::io::Read{
    let mut lock = reader.lock().unwrap();
    lock.next()
}

pub fn map<S>(seqs: DnaString,
              dbg: &DebruijnGraph<KmerType, EqClassIdType>,
              eq_classes: &Vec<Vec<S>>) -> Vec<S>
where S: Ord + PartialEq + Clone {
    let mut eq_class: Vec<S> = Vec::new();

    for kmer in seqs.iter_kmers() {
        let (nid, _, _) = match dbg.find_link(kmer, Dir::Right){
            Some(links) => links,
            None => (std::usize::MAX, Dir::Right, false),
        };

        if nid != std::usize::MAX {
            let eq_id = *dbg.get_node(nid).data();
            let labels = &eq_classes[eq_id as usize];
            eq_class.extend(labels.clone().into_iter());
            pdqsort::sort(&mut eq_class);
            eq_class.dedup();
        }
        //let maybe_data = phf.get(&kmer);
        //match maybe_data {
        //    Some((_, color)) => {
        //        let labels = eq_classes[*color as usize].clone();
        //        eq_class.extend(labels) ;
        //        pdqsort::sort(&mut eq_class);
        //        eq_class.dedup();
        //    },
        //    None => (),
        //}
        //println!("{:?}", eq_class);
    }

    eq_class
}
