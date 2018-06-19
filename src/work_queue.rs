// Part of code taken from
// https://gist.github.com/LeoTindall/e6d40782b05dc8ac40faf3a0405debd3
use std;
use std::hash::Hash;
use std::fmt::Debug;
use std::sync::{Mutex, Arc};
use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};


use boomphf;
//use pdqsort;
use debruijn::*;
use debruijn::graph::*;
use debruijn::filter::*;
use bio::io::{fastq};
use debruijn::compression::*;
use docks::generate_msps;
use debruijn::dna_string::{DnaString, DnaStringSlice};
use config::{KmerType, STRANDED, REPORT_ALL_KMER, MEM_SIZE, L, DocksUhs};

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
               -> (std::vec::Vec<(u16, DnaStringSlice<'a>, Exts)>, usize){
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
                if bucket_id >= uhs.len() as u16{
                    panic!("Bucket size: {:?} id: {:?} out of bound",
                           uhs.len(), bucket_id);
                }

                let slice = seq.slice(msp.start(), msp.end());
                let exts = Exts::from_dna_string(seq, msp.start(), msp.len());
                bucket_slices.push((bucket_id, slice, exts));
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

pub fn map<S>(read_seq: DnaString,
              dbg: &DebruijnGraph<KmerType, EqClassIdType>,
              eq_classes: &Vec<Vec<S>>) -> Vec<S>
where S: Clone + Ord + PartialEq + Debug + Sync + Send + Hash {
    let mut all_colors: Vec<Vec<S>> = Vec::new();
    let read_length = read_seq.len();

    let mut kmer_pos = 0;
    let mut kmer = read_seq.first_kmer();
    let kmer_length = KmerType::k();
    let last_kmer_pos = read_length - kmer_length;

    while kmer_pos <= last_kmer_pos {
        let nid = match dbg.find_link(kmer, Dir::Right){
            Some((nid, _, _)) => nid,
            None => {
                //match dbg.find_link(kmer, Dir::Left) {
                //    Some(links) => links,
                //    None => (std::usize::MAX, Dir::Left, true),
                //}
                std::usize::MAX
            },
        };

        if nid != std::usize::MAX {
            // get the node
            let ref_node = dbg.get_node(nid);
            let ref_seq_slice = ref_node.sequence();
            let ref_length = ref_seq_slice.len();

            // useless sanity check but keeping it
            let ref_kmer = ref_seq_slice.first_kmer();
            if kmer != ref_kmer {
                panic!("{:?} should be equal to {:?}",
                       kmer, ref_kmer);
            }

            let read_offset = kmer_pos + kmer_length;
            let remaining_read = read_length - read_offset;
            // TODO: verify if it's always going to match at the start
            // of the node
            let informative_ref = ref_length - kmer_length;
            let max_matchable_pos = std::cmp::min(remaining_read, informative_ref);

            for idx in 0..max_matchable_pos {
                let ref_pos = kmer_length + idx;
                let read_pos = read_offset + idx;

                // compare base by base
                if ref_seq_slice.get(ref_pos) != read_seq.get(read_pos) {
                    break;
                }
                kmer_pos += 1;
            }

            // extract colors
            let eq_id = ref_node.data();
            let colors = &eq_classes[*eq_id as usize];

            all_colors.push(colors.clone());

            println!("{:?}, {:?}, {:?}", ref_seq_slice, ref_node, colors);
        } // end-if

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

        kmer_pos += 1;
        if kmer_pos > last_kmer_pos {
            break;
        }
        kmer = read_seq.get_kmer(kmer_pos);
    }// end-loop

    // Take the intersection of the sets
    let total_classes = all_colors.len();
    if total_classes == 0 {
        Vec::new()
    }
    else{
        let elem: Vec<S> = all_colors.pop().unwrap();
        if total_classes == 1 {
            return elem
        }

        let mut eq_class_set: HashSet<S> = elem.into_iter().collect();
        for colors in all_colors {
            let colors_set: HashSet<S> = colors.into_iter().collect();
            eq_class_set = eq_class_set.intersection(&colors_set).cloned().collect();
        }

        let eq_class: Vec<S> = eq_class_set.into_iter().collect();
        eq_class
    }
}
