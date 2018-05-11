extern crate debruijn;

use std::hash::Hash;
use std::marker::PhantomData;
use std::collections::HashSet;

use debruijn::Exts;
use debruijn::{filter, compression};

// Extending trait KmerSummarizer for fitering
pub struct ScmapCountFilterSet<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>,
}

impl<D> ScmapCountFilterSet<D> {
    pub fn new(min_kmer_obs: usize) -> ScmapCountFilterSet<D> {
        ScmapCountFilterSet {
            min_kmer_obs: min_kmer_obs,
            phantom: PhantomData,
        }
    }
}

impl<D: Ord + Hash> filter::KmerSummarizer<D, HashSet<D>> for ScmapCountFilterSet<D> {
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, HashSet<D>) {
        let mut nobs = 0;
        let mut all_exts = Exts::empty();
        let mut out_data: HashSet<D> = HashSet::new();

        for (_, exts, d) in items {
            out_data.insert(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        (nobs as usize >= self.min_kmer_obs, all_exts, out_data)
    }
}


// Extending trait CompressionSpec for compression
pub struct ScmapCompress<D, F> {
    func: F,
    d: PhantomData<D>,
}

impl<D, F> ScmapCompress<D, F> {
    pub fn new(func: F) -> ScmapCompress<D, F> {
        ScmapCompress {
            func: func,
            d: PhantomData,
        }
    }
}

impl<D: PartialEq, F> compression::CompressionSpec<D> for ScmapCompress<D, F>
where
    for<'r> F: Fn(D, &'r D) -> D,
{
    fn reduce(&self, d: D, other: &D) -> D {
        (self.func)(d, other)
    }

    fn join_test(&self, d1: &D, d2: &D) -> bool {
        if d1 == d2 { true } else { false }
    }
}
