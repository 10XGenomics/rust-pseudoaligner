// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use std;
use std::cmp::Ordering;

use boomphf::hashmap::{NoKeyBoomHashMap};
use debruijn::{Dir, Mer, Kmer, Vmer};
use debruijn::dna_string::{DnaString};
use debruijn::graph::DebruijnGraph;
use debruijn::filter::EqClassIdType;


#[derive(Serialize, Deserialize, Debug)]
pub struct Pseudoaligner<K: Kmer>  
{
    dbg: DebruijnGraph<K, EqClassIdType>,
    eq_classes: Vec<Vec<u32>>,
    dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
}

impl<K> Pseudoaligner<K>
    where K: Kmer + Sync + Send
{
    pub fn new(dbg: DebruijnGraph<K, EqClassIdType>, 
            eq_classes: Vec<Vec<u32>>,
            dbg_index: NoKeyBoomHashMap<K, (u32, u32)>) -> Pseudoaligner<K>
    {
        Pseudoaligner { dbg, eq_classes, dbg_index }
    }
    /// Pseudo-align `read_seq` to determine its the equivalence class.
    pub fn map_read(&self, read_seq: &DnaString) -> Option<(Vec<u32>, usize)>
    {
        let read_length = read_seq.len();
        let mut read_coverage: usize = 0;
        let mut colors: Vec<u32> = Vec::new();
        let left_extend_threshold = (0.4 * read_length as f32) as usize;

        let mut kmer_pos: usize = 0;
        let kmer_length = K::k();
        let last_kmer_pos = read_length - kmer_length;

        // Scan the read for the first kmer that exists in the reference
        let find_kmer_match = |kmer_pos: &mut usize| -> Option<(usize, usize)> {
            while *kmer_pos <= last_kmer_pos {
                let read_kmer = read_seq.get_kmer(*kmer_pos);

                match self.dbg_index.get(&read_kmer) {
                    None => (),
                    Some((nid, offset)) => {
                        let node = self.dbg.get_node(*nid as usize);
                        let ref_seq_slice = node.sequence();
                        let ref_kmer: K = ref_seq_slice.get_kmer(*offset as usize);

                        if read_kmer == ref_kmer {
                            return Some((*nid as usize, *offset as usize));
                        }
                    }
                };
                *kmer_pos += 1;
            }

            None
        };

        // extract the first exact matching position of read
        let (mut node_id, mut kmer_offset) = 
            // get the first match through mphf
            match find_kmer_match(&mut kmer_pos) {
                None => (None, None),
                Some((nid, offset)) => (Some(nid), Some(offset))
            };

        // check if we can extend back if there were SNP in every kmer query
        if kmer_pos >= left_extend_threshold && node_id.is_some() {
            let mut last_pos = kmer_pos - 1;
            let mut prev_node_id = node_id.unwrap();
            let mut prev_kmer_offset = kmer_offset.unwrap() - 1;

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
                    let color = prev_node.data();
                    colors.push(*color);
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
                    if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
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

        // Take the intersection of the sets
        let colors_len = colors.len();
        if colors_len == 0 {
            if read_coverage != 0 {
                panic!(
                    "Different read coverage {:?} than num of eqclasses {:?}",
                    colors_len, read_coverage
                );
            }

            None
        } else {

            // Intersect the equivalence classes
            let first_color = colors.pop().unwrap();
            let mut eq_class = self.eq_classes[first_color as usize].clone();

            for color in colors {
                intersect(&mut eq_class, &self.eq_classes[color as usize]);
            }

            Some((eq_class, read_coverage))
        }
    }
}

/// Compute the intersection of v1 and v2 inplace on top of v1
/// v1 and v2 must be sorted
fn intersect<T: Eq + Ord>(v1: &mut Vec<T>, v2: &[T]) {
    if v1.len() == 0 {
        return;
    }

    if v2.len() == 0 {
        v1.clear();
    }

    let mut fill_idx1 = 0;
    let mut idx1 = 0;
    let mut idx2 = 0;

    while idx1 < v1.len() && idx2 < v2.len() {
        match v1[idx1].cmp(&v2[idx2]) {
            Ordering::Less => idx1 += 1,
            Ordering::Greater => idx2 += 1,
            Ordering::Equal => { 
                v1.swap(fill_idx1, idx1);
                idx1 += 1;
                idx2 += 1;
                fill_idx1 += 1;
            }
        }
    }
}

