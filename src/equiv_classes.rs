// Copyright 2017 10x Genomics

//! Generate equivalence classes for pseudoaligner
use std::fmt::Debug;
use std::hash::Hash;
use std::ops::Deref;
use std::sync::atomic::{AtomicUsize, Ordering};

use dashmap::DashMap;

use debruijn::filter::KmerSummarizer;
use debruijn::Exts;

//Equivalence class based implementation
pub type EqClassIdType = u32;
pub struct CountFilterEqClass<D: Eq + Hash + Send + Sync + Debug + Clone> {
    min_kmer_obs: usize,
    eq_classes: DashMap<Vec<D>, EqClassIdType>,
    num_eq_classes: AtomicUsize,
}

impl<D: Eq + Hash + Send + Sync + Debug + Clone> CountFilterEqClass<D> {
    pub fn new(min_kmer_obs: usize) -> CountFilterEqClass<D> {
        CountFilterEqClass {
            min_kmer_obs: min_kmer_obs,
            eq_classes: DashMap::<Vec<D>, EqClassIdType>::new(4),
            num_eq_classes: AtomicUsize::new(0),
        }
    }

    pub fn get_eq_classes(&self) -> Vec<Vec<D>> {
        let mut eq_class_vec = Vec::new();
        eq_class_vec.resize(self.get_number_of_eq_classes(), Vec::new());

        let mut eq_ids = Vec::new();

        for item in self.eq_classes.iter() {
            eq_class_vec[*item.value() as usize] = item.key().clone();
            eq_ids.push(*item.value() as usize)
        }

        // consistency property the equivalence classes must be assigned
        // unique ids from 0 to N, with no gaps. This could be violated
        // if theres is a race condition when assigning equivalence class
        // ids in CountFilterEqClass::summarize below.  panic if this
        // property doesn't hold.
        eq_ids.sort();
        for i in 0..eq_ids.len() {
            assert_eq!(eq_ids[i], i);
        }

        eq_class_vec
    }

    pub fn get_number_of_eq_classes(&self) -> usize {
        self.num_eq_classes.load(Ordering::SeqCst)
    }

    pub fn fetch_add(&self) -> usize {
        self.num_eq_classes.fetch_add(1, Ordering::SeqCst)
    }
}

impl<D: Eq + Ord + Hash + Send + Sync + Debug + Clone> KmerSummarizer<D, EqClassIdType>
    for CountFilterEqClass<D>
{
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(
        &self,
        items: F,
    ) -> (bool, Exts, EqClassIdType) {
        let mut all_exts = Exts::empty();

        // the ids of the sequences in the equivalence class
        let mut eq_class = Vec::new();

        let mut nobs = 0;
        for (_, exts, d) in items {
            eq_class.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        eq_class.sort();
        eq_class.dedup();

        // register the equivalence class and assign it a unique id.
        // IDs must be sequential from 0 to N. This must be an atomic operation.
        // The correctness of the eqclass_ids is checked above int get_eq_classes.
        let eq_ref = self.eq_classes.get_or_insert_with(&eq_class, || {
            self.num_eq_classes.fetch_add(1, Ordering::SeqCst) as u32
        });

        let eq_id = eq_ref.deref().clone();
        (nobs as usize >= self.min_kmer_obs, all_exts, eq_id)
    }
}
