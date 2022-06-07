// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use std::collections::HashMap;
use std::fmt::Debug;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{mpsc, Arc, Mutex};
use std::{self, fs::File, str};

use bio::io::fastq;
use boomphf::hashmap::NoKeyBoomHashMap;
use crossbeam_utils::thread::scope;
use debruijn::dna_string::DnaString;

use anyhow::Error;
use debruijn::graph::DebruijnGraph;
use debruijn::{Dir, Kmer, Mer, Vmer};
use log::info;
use serde::{Deserialize, Serialize};

use crate::config::{DEFAULT_ALLOWED_MISMATCHES, LEFT_EXTEND_FRACTION, READ_COVERAGE_THRESHOLD};
use crate::equiv_classes::EqClassIdType;
use crate::utils;

#[derive(Serialize, Deserialize, Debug)]
pub struct Pseudoaligner<K: Kmer> {
    pub dbg: DebruijnGraph<K, EqClassIdType>,
    pub eq_classes: Vec<Vec<u32>>,
    pub dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
    pub tx_names: Vec<String>,
    pub tx_gene_mapping: HashMap<String, String>,
}

impl<K: Kmer + Sync + Send> Pseudoaligner<K> {
    pub fn new(
        dbg: DebruijnGraph<K, EqClassIdType>,
        eq_classes: Vec<Vec<u32>>,
        dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
        tx_names: Vec<String>,
        tx_gene_mapping: HashMap<String, String>,
    ) -> Pseudoaligner<K> {
        Pseudoaligner {
            dbg,
            eq_classes,
            dbg_index,
            tx_names,
            tx_gene_mapping,
        }
    }

    /// Pseudo-align `read_seq` and return a list of nodes that the read was aligned to, with mismatch = 2
    pub fn map_read_to_nodes(&self, read_seq: &DnaString, nodes: &mut Vec<usize>) -> Option<usize> {
        self.map_read_to_nodes_with_mismatch(read_seq, nodes, DEFAULT_ALLOWED_MISMATCHES).map(|(read_coverage, _mismatches)| read_coverage)
    }

    /// Pseudo-align `read_seq` and return a list of nodes that the read was aligned to, with configurable # of allowed mismatches
    pub fn map_read_to_nodes_with_mismatch(
        &self,
        read_seq: &DnaString,
        nodes: &mut Vec<usize>,
        allowed_mismatches: usize,
    ) -> Option<(usize, usize)> {
        let read_length = read_seq.len();
        let mut read_coverage: usize = 0;
        let mut mismatch_count: usize = 0;

        // We're filling out nodes
        nodes.clear();

        let left_extend_threshold = (LEFT_EXTEND_FRACTION * read_length as f64) as usize;

        let mut kmer_pos: usize = 0;
        let kmer_length = K::k();

        if read_seq.len() < kmer_length {
            return None;
        }

        let last_kmer_pos = read_length - kmer_length;
        let mut kmer_lookups = 0;

        {
            // Scan the read for the first kmer that exists in the reference
            let mut find_kmer_match = |kmer_pos: &mut usize| -> Option<(usize, usize)> {
                while *kmer_pos <= last_kmer_pos {
                    let read_kmer = read_seq.get_kmer(*kmer_pos);

                    kmer_lookups += 1;
                    match self.dbg_index.get(&read_kmer) {
                        None => (),
                        Some((nid, offset)) => {
                            // Verify that the kmer actually matches -- the MPHF can have false
                            // positives.
                            let node = self.dbg.get_node(*nid as usize);
                            let ref_seq_slice = node.sequence();
                            let ref_kmer: K = ref_seq_slice.get_kmer(*offset as usize);

                            if read_kmer == ref_kmer {
                                return Some((*nid as usize, *offset as usize));
                            }
                        }
                    };
                    *kmer_pos += 3;
                }

                None
            };

            // extract the first exact matching position of a kmer
            // from the read in the DBG
            let (mut node_id, mut kmer_offset) = match find_kmer_match(&mut kmer_pos) {
                None => (None, None),
                Some((nid, offset)) => (Some(nid), Some(offset)),
            };

            // check if we can extend back if there were SNP in every kmer query
            if kmer_pos >= left_extend_threshold && node_id.is_some() {
                let mut last_pos = kmer_pos - 1;
                let mut prev_node_id = node_id.unwrap();
                let mut prev_kmer_offset = if kmer_offset.unwrap() > 0 {
                    kmer_offset.unwrap() - 1
                } else {
                    0
                };

                loop {
                    let node = self.dbg.get_node(prev_node_id);
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
                        if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
                            // Record mismatch
                            mismatch_count += 1;

                            // Allowing num_mismatch-SNP
                            seen_snp += 1;
                            if seen_snp > allowed_mismatches {
                                premature_break = true;
                                break;
                            }
                        }

                        matched_bases += 1;
                        read_coverage += 1;
                    }

                    //break the loop if end of read reached or a premature mismatch
                    if last_pos + 1 - matched_bases == 0 || premature_break {
                        break;
                    }

                    // adjust last position
                    last_pos -= matched_bases;

                    // If reached here then a fork is found in the reference.
                    let exts = node.exts();
                    let next_base = read_seq.get(last_pos);
                    if exts.has_ext(Dir::Left, next_base) {
                        // found a left extention.
                        let index = exts
                            .get(Dir::Left)
                            .iter()
                            .position(|&x| x == next_base)
                            .unwrap();

                        let edge = node.l_edges()[index];

                        //update the previous node's id
                        prev_node_id = edge.0;
                        let prev_node = self.dbg.get_node(prev_node_id);
                        prev_kmer_offset = prev_node.sequence().len() - kmer_length;

                        // extract colors
                        nodes.push(prev_node.node_id);
                    } else {
                        break;
                    }
                } // end-loop
            } //end-if

            // forward search
            if kmer_pos <= last_kmer_pos {
                loop {
                    let node = self.dbg.get_node(node_id.unwrap());
                    //println!("{:?}, {:?}, {:?}, {:?}",
                    //         node, node.sequence(),
                    //         &eq_classes[ *node.data() as usize],
                    //         kmer_offset);
                    kmer_pos += kmer_length;
                    read_coverage += kmer_length;

                    // extract colors
                    nodes.push(node.node_id);

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
                        if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
                            // Record mismatch
                            mismatch_count += 1;

                            // Allowing num_mismatch-SNP
                            seen_snp += 1;
                            if seen_snp > allowed_mismatches {
                                premature_break = true;
                                break;
                            }
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
                    let next_base = read_seq.get(kmer_pos);

                    if !premature_break && exts.has_ext(Dir::Right, next_base) {
                        // found a right extention.
                        let index = exts
                            .get(Dir::Right)
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
                    } else {
                        // can't extend node in dbg extract read using mphf
                        // TODO: might have to check some cases
                        if kmer_pos > last_kmer_pos {
                            // can't search in mphf if no full kmer can be made
                            break;
                        }

                        // get the match through mphf
                        match find_kmer_match(&mut kmer_pos) {
                            None => break,
                            Some((nid, offset)) => {
                                node_id = Some(nid);
                                kmer_offset = Some(offset);
                            }
                        };
                    }
                } // end-loop
            } //end-if
        }

        if nodes.is_empty() {
            if read_coverage != 0 {
                panic!(
                    "Different read coverage {:?} than num of eqclasses {:?}",
                    nodes.len(),
                    read_coverage
                );
            }
            //println!("lookups: {} -- no hit", kmer_lookups);
            None
        } else {
            //println!("lookups: {} -- cov: {}", kmer_lookups, read_coverage);
            Some((read_coverage, mismatch_count))
        }
    }

    /// Convert a list of nodes contacted by a read into an equivalence class.
    /// Supply node list in `nodes`. Equivalence class will be written to `eq_class`.
    pub fn nodes_to_eq_class(&self, nodes: &mut Vec<usize>, eq_class: &mut Vec<u32>) {
        eq_class.clear();

        if nodes.is_empty() {
            return;
        }

        // Sort nodes to get the shorter equivalence class first.
        nodes.sort_by_key(|n| {
            let eqclass_id = self.dbg.get_node(*n).data();
            self.eq_classes[*eqclass_id as usize].len()
        });

        let _lens: Vec<_> = nodes
            .iter()
            .map(|n| {
                let eqclass_id = self.dbg.get_node(*n).data();
                self.eq_classes[*eqclass_id as usize].len()
            })
            .collect();
        //println!("nodes: {:?}, lens: {:?}", nodes, lens);

        // Intersect the equivalence classes
        let first_node = nodes[0];

        //println!("node: {}, seq: {:?}", first_node,  self.dbg.get_node(first_node).sequence());
        let first_color = self.dbg.get_node(first_node).data();
        eq_class.extend(&self.eq_classes[*first_color as usize]);

        for node in nodes.iter().skip(1) {
            let color = self.dbg.get_node(*node).data();
            intersect(eq_class, &self.eq_classes[*color as usize]);
        }
    }

    /// Pseudoalign the `read_seq` to the graph. Returns a tuple of the
    /// eqivalence class, the number of bases aligned on success,
    /// and the number of mismatched bases, or None is no alignment could be found.
    pub fn map_read_with_mismatch(
        &self,
        read_seq: &DnaString,
        allowed_mismatches: usize,
    ) -> Option<(Vec<u32>, usize, usize)> {
        let mut nodes = Vec::new();

        match self.map_read_to_nodes_with_mismatch(read_seq, &mut nodes, allowed_mismatches) {
            Some((read_coverage, mismatches)) => {
                let mut eq_class = Vec::new();
                self.nodes_to_eq_class(&mut nodes, &mut eq_class);
                Some((eq_class, read_coverage, mismatches))
            }
            None => None,
        }
    }

    /// Pseudoalign the `read_seq` to the graph with # mismatches = 2. Returns a tuple of the
    /// eqivalence class and the number of bases aligned on success
    /// or None is no alignment could be found.
    pub fn map_read(&self, read_seq: &DnaString) -> Option<(Vec<u32>, usize)> {
        self.map_read_with_mismatch(read_seq, DEFAULT_ALLOWED_MISMATCHES).map(|(eq_class, read_coverage, _mismatches)| (eq_class, read_coverage))
    }
}

/// Compute the intersection of v1 and v2 inplace on top of v1
/// v1 and v2 must be sorted and deduplicated.
pub fn intersect<T: Eq + Ord>(v1: &mut Vec<T>, v2: &[T]) {
    if v1.is_empty() {
        return;
    }

    if v2.is_empty() {
        v1.clear();
    }

    let mut fill_idx1 = 0;
    let mut idx1 = 0;
    let mut idx2 = 0;

    while idx1 < v1.len() && idx2 < v2.len() {
        let rem_slice = &v2[idx2..];
        match rem_slice.binary_search(&v1[idx1]) {
            Ok(pos) => {
                v1.swap(fill_idx1, idx1);
                fill_idx1 += 1;
                idx1 += 1;
                idx2 = pos + 1;
            }
            Err(pos) => {
                idx1 += 1;
                idx2 = pos;
            }
        }
    }
    v1.truncate(fill_idx1);
}

pub fn process_reads<K: Kmer + Sync + Send, P: AsRef<Path> + Debug>(
    reader: fastq::Reader<io::BufReader<File>>,
    index: &Pseudoaligner<K>,
    outdir: P,
    num_threads: usize,
) -> Result<(), Error> {
    info!("Done Reading index");
    info!("Starting Multi-threaded Mapping");
    info!("Output directory: {:?}", outdir);

    let (tx, rx) = mpsc::sync_channel(num_threads);
    let atomic_reader = Arc::new(Mutex::new(reader.records()));

    info!("Spawning {} threads for Mapping.\n", num_threads);
    scope(|scope| {
        for _ in 0..num_threads {
            let tx = tx.clone();
            let reader = Arc::clone(&atomic_reader);

            scope.spawn(move |_| {
                loop {
                    // If work is available, do that work.
                    match utils::get_next_record(&reader) {
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
                    if dead_thread_count == num_threads {
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
    })
    .unwrap(); //end crossbeam

    eprintln!();
    info!("Done Mapping Reads");
    Ok(())
}

#[cfg(test)]
mod test {
    use super::*;
    use proptest::collection::vec;
    use proptest::prelude::*;
    use proptest::proptest;
    use std::collections::HashSet;
    use std::hash::Hash;
    use std::iter::FromIterator;

    fn test_intersect<T: Hash + Eq + Clone + Ord + Debug>(v1: &Vec<T>, v2: &Vec<T>) {
        let mut c1 = v1.clone();
        let c2 = v2.clone();

        let s1: HashSet<T> = HashSet::from_iter(c1.iter().cloned());
        let s2: HashSet<T> = HashSet::from_iter(c2.iter().cloned());
        let intersection = s1.intersection(&s2);

        let mut int1: Vec<T> = intersection.cloned().collect();
        int1.sort();

        intersect(&mut c1, &c2);

        assert_eq!(c1, int1);
    }

    #[test]
    fn intersect_test() {
        let v1 = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let v2 = vec![1, 2, 3];
        let v3 = vec![1, 4, 5];
        let v4 = vec![7, 8, 9];
        let v5 = vec![9];
        let v6: Vec<usize> = vec![];
        let v7 = vec![1, 2, 3, 6, 7, 8, 9];
        let v8 = vec![1, 7, 8, 9, 10];
        let v9 = vec![10, 15, 20];
        let v10 = vec![21, 22, 23];
        let v11 = vec![0];
        let v12 = vec![0, 1000, 5000];
        let v13 = vec![0, 1000, 1000001];
        let v14 = vec![5];
        let v15 = vec![100000000];
        let v16 = vec![1, 23, 45, 1000001, 100000000];

        let vecs = vec![
            v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16,
        ];

        for v1 in vecs.iter() {
            for v2 in vecs.iter() {
                test_intersect(v1, v2);
                test_intersect(v2, v1);
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 1000, .. ProptestConfig::default()})]
        #[test]
        fn intersect_prop_test(
            mut v1 in vec(0..100usize, 0..5000usize),
            mut v2 in vec(0..100usize, 0..5000usize),
        ) {

            v1.sort_unstable(); v1.dedup();
            v2.sort_unstable(); v2.dedup();
            test_intersect(&v1, &v2);
            test_intersect(&v2, &v1);
        }
    }
}
