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

pub fn analyze<S>( bucket_data: Vec<(DnaStringSlice, Exts, S)>,
                   summarizer: &Arc<CountFilterEqClass<S>>)
                   -> Option<BaseGraph<KmerType, EqClassIdType>>
where S: Clone + Eq + Hash + Ord + Debug + Send + Sync {
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
    let combined_graph = BaseGraph::combine(uncompressed_dbgs.into_iter())
        .finish();

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
              eq_classes: &Vec<Vec<S>>,
              phf: &boomphf::NoKeyBoomHashMap2<KmerType, usize, u32>)
              -> Option<(Vec<S>, usize)>
where S: Clone + Ord + PartialEq + Debug + Sync + Send + Hash {
    let read_length = read_seq.len();
    let mut read_coverage: usize = 0;
    let mut colors: Vec<u32> = Vec::new();
    let left_extend_threshold = (0.4 * read_length as f32 ) as usize;

    let mut kmer_pos: usize = 0;
    let kmer_length = KmerType::k();
    let last_kmer_pos = read_length - kmer_length;

    // extract the first exact matching position of read
    let mut node_id = None;
    let mut kmer_offset = None;

    let get_data = | kmer_pos: &mut usize |
                                  -> Option<(usize, usize)> {
        while *kmer_pos <= last_kmer_pos {
            let read_kmer = read_seq.get_kmer(*kmer_pos);

            match phf.get(&read_kmer) {
                None => (),
                Some((&nid, &offset)) => {

                    let node = dbg.get_node( nid );
                    let ref_seq_slice = node.sequence();

                    let read_kmer = read_seq.get_kmer::<KmerType>(*kmer_pos);
                    let ref_kmer = ref_seq_slice.get_kmer::<KmerType>(offset as usize);

                    if read_kmer == ref_kmer {
                        return Some((nid, offset as usize));
                    }
                }
            };
            *kmer_pos += 1;
        }

        return None;
    };

    // get the first match through mphf
    match get_data(&mut kmer_pos) {
        None => (),
        Some((nid, offset)) => {
            node_id = Some(nid);
            kmer_offset = Some(offset);
        },
    };

    // check if we can extend back if there were SNP in every kmer query
    if kmer_pos >= left_extend_threshold && ! node_id.is_none() {
        let mut last_pos = kmer_pos - 1;
        let mut prev_node_id = node_id.unwrap();
        let mut prev_kmer_offset = kmer_offset.unwrap() - 1;

        loop {
            let node = dbg.get_node(prev_node_id);
            //println!("{:?}, {:?}, {:?}, {:?}, {:?}",
            //         node, node.sequence(),
            //         &eq_classes[ *node.data() as usize],
            //         prev_kmer_offset, last_pos);


            // length of remaining read before kmer match
            let skipped_read = last_pos + 1;

            // length of the skipped node sequence before kmer match
            let skipped_ref = prev_kmer_offset + 1;

            // find maximum extention possbile before fork or eof read
            let max_matchable_pos = std::cmp::min(skipped_read, skipped_ref);

            let ref_seq_slice = node.sequence();
            let mut premature_break = false;
            let mut matched_bases = 0;
            let mut seen_snp = 0;
            for idx in 0..max_matchable_pos {
                let ref_pos = prev_kmer_offset - idx;
                let read_offset = last_pos - idx;

                // compare base by base
                if ref_seq_slice.get( ref_pos ) != read_seq.get( read_offset ) {
                    if seen_snp > 3 {
                        premature_break = true;
                        break;
                    }

                    // Allowing 2-SNP
                    seen_snp += 1;
                }

                matched_bases += 1;
                read_coverage += 1;
            }


            //break the loop if end of read reached or a premature mismatch
            if last_pos - matched_bases + 1 == 0 || premature_break {
                break;
            }

            // adjust last position
            last_pos -= matched_bases;

            // If reached here then a fork is found in the reference.
            let exts = node.exts();
            let next_base = read_seq.get( last_pos );
            if exts.has_ext(Dir::Left, next_base) {
                // found a left extention.
                let index = exts.get(Dir::Left)
                    .iter()
                    .position(|&x| x == next_base)
                    .unwrap();

                let edge = node.l_edges()[index];

                //update the previous node's id
                prev_node_id = edge.0;
                let prev_node = dbg.get_node(prev_node_id);
                prev_kmer_offset = prev_node.sequence().len() - kmer_length;

                // extract colors
                let color = prev_node.data();
                colors.push(*color);
            }
            else {
                break;
            }
        } // end-loop
    } //end-if

    // forward search
    if kmer_pos <= last_kmer_pos {
        loop {
            let node = dbg.get_node(node_id.unwrap());
            //println!("{:?}, {:?}, {:?}, {:?}",
            //         node, node.sequence(),
            //         &eq_classes[ *node.data() as usize],
            //         kmer_offset);
            kmer_pos += kmer_length;
            read_coverage += kmer_length;

            // extract colors
            let color = node.data();
            colors.push(*color);

            // length of remaining read after kmer match
            let remaining_read = read_length - kmer_pos;

            // length of the remaining node sequence after kmer match
            let ref_seq_slice = node.sequence();
            let ref_length = ref_seq_slice.len();
            let ref_offset = kmer_offset.unwrap() + kmer_length;
            let informative_ref = ref_length - ref_offset;

            // find maximum extention possbile before fork or eof read
            let max_matchable_pos = std::cmp::min(remaining_read, informative_ref);

            let mut premature_break = false;
            let mut matched_bases = 0;
            let mut seen_snp = 0;
            for idx in 0..max_matchable_pos {
                let ref_pos = ref_offset + idx;
                let read_offset = kmer_pos + idx;

                // compare base by base
                if ref_seq_slice.get( ref_pos ) != read_seq.get( read_offset ) {
                    if seen_snp > 3 {
                        premature_break = true;
                        break;
                    }

                    // Allowing 2-SNP
                    seen_snp += 1;
                }

                matched_bases += 1;
                read_coverage += 1;
            }

            kmer_pos += matched_bases;
            //break the loop if end of read reached or a premature mismatch
            if kmer_pos >= read_length {
                break;
            }

            // If reached here then a fork is found in the reference.
            let exts = node.exts();
            let next_base = read_seq.get( kmer_pos );

            if !premature_break && exts.has_ext(Dir::Right, next_base) {
                // found a right extention.
                let index = exts.get(Dir::Right)
                    .iter()
                    .position(|&x| x == next_base)
                    .unwrap();

                let edge = node.r_edges()[index];

                //update the next node's id
                node_id = Some(edge.0);
                kmer_offset = Some(0);

                //adjust for kmer_position
                kmer_pos -= kmer_length - 1;
                read_coverage -= kmer_length - 1;
            }
            else{
                // can't extend node in dbg extract read using mphf
                // TODO: might have to check some cases
                if kmer_pos > last_kmer_pos {
                    // can't search in mphf if no full kmer can be made
                    break;
                }

                // get the match through mphf
                match get_data(&mut kmer_pos) {
                    None => break,
                    Some((nid, offset)) => {
                        node_id = Some(nid);
                        kmer_offset = Some(offset);
                    },
                };
            }
        } // end-loop
    }//end-if

    // Take the intersection of the sets
    let colors_len = colors.len();
    if colors_len == 0 {
        if read_coverage != 0 {
            panic!("Different read coverage {:?} than num of eqclasses {:?}",
                   colors_len, read_coverage);
        }

        return None
    }
    else{
        let color: u32 = colors.pop().unwrap();
        let eq_class: Vec<S> = eq_classes[color as usize].to_owned();

        if colors_len == 1 {
            return Some((eq_class, read_coverage))
        }

        let mut eq_class_set: HashSet<S> = eq_class.into_iter().collect();
        for color in colors {
            let eq_class: HashSet<S> = eq_classes[color as usize]
                .to_owned()
                .into_iter()
                .collect();

            eq_class_set = eq_class_set
                .intersection(&eq_class)
                .cloned()
                .collect();
        }

        let eq_classes_vec: Vec<S> = eq_class_set.into_iter().collect();
        return Some((eq_classes_vec, read_coverage))
    }
}
