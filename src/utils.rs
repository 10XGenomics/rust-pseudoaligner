use std;
use std::mem;
use std::hash::Hash;
use std::fmt::Debug;
use std::marker::PhantomData;
use std::collections::{HashMap, HashSet};

use heapsize::HeapSizeOf;
use pdqsort;
use boomphf;
use std::str::FromStr;
use bit_set::BitSet;
use itertools::Itertools;
use debruijn::{compression, filter};
use debruijn::{Dir, Kmer, Exts, Vmer, Mer};
use bio::data_structures::smallints::SmallInts;

// Extending trait CompressionSpec for compression
pub struct ScmapCompress<D> {
    d: PhantomData<D>,
}

impl<D> ScmapCompress<D> {
    pub fn new() -> ScmapCompress<D> {
        ScmapCompress {
            d: PhantomData,
        }
    }
}

impl<D: PartialEq> compression::CompressionSpec<D> for ScmapCompress<D>
where D: std::fmt::Debug
{
    fn reduce(&self, d: D, other: &D) -> D {
        if d != *other {
            panic!("{:?} != {:?}, Should not happen", d, *other);
        }
        d
    }

    fn join_test(&self, d1: &D, d2: &D) -> bool {
        if d1 == d2 { true } else { false }
    }
}

pub struct CountFilterSmallInt<D> {
    min_kmer_obs: usize,
    phantom: PhantomData<D>,
}

impl<D> CountFilterSmallInt<D> {
    pub fn new(min_kmer_obs: usize) -> CountFilterSmallInt<D> {
        CountFilterSmallInt {
            min_kmer_obs: min_kmer_obs,
            phantom: PhantomData,
        }
    }
}


impl filter::KmerSummarizer<usize, SmallInts<u8, usize>> for CountFilterSmallInt<usize> {
    fn summarize<K, F: Iterator<Item = (K, Exts, usize)>>(&self, items: F) -> (bool, Exts, SmallInts<u8, usize>) {
        let mut all_exts = Exts::empty();
        let mut out_data: Vec<usize> = Vec::new();

        let mut nobs = 0;
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();  out_data.dedup();

        let mut deduplicated_data: SmallInts<u8, usize> = SmallInts::new();
        for value in out_data{
            deduplicated_data.push(value);
        }

        (nobs as usize >= self.min_kmer_obs, all_exts, deduplicated_data)
    }
}


// Minimal perfect hash
pub struct BoomHashMap<K: Clone + Hash + Debug, Exts, D> {
    mph_hash: boomphf::Mphf<K>,
    keys: Vec<K>,
    exts: Vec<Exts>,
    data: Vec<D>
}

/// Read a shard and determine the valid kmers
/// Low memory implementation that should consume < 4G of temporary memory
/// To reduce memory consumption, set track_bcs to false to forget about BC lists.
#[inline(never)]
pub fn filter_kmers_with_mphf<K: Kmer, V: Vmer<K>, D1: Clone, DS, S: filter::KmerSummarizer<D1, DS>>(
    seqs: Vec<V>,
    summarizer: S,
    stranded: bool,
    report_all_kmers: bool,
    memory_size: usize,
){// -> BoomHashMap<K, Exts, DS> {

    let rc_norm = !stranded;

    // let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    // let mut valid_data = Vec::new();

    let num_seqs: usize = seqs.len();

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs.iter()
        .map(|ref vmer| vmer.len().saturating_sub(K::k() - 1))
        .sum();
    let kmer_mem = input_kmers * mem::size_of::<(K, D1)>();
    let max_mem = memory_size * (10 as usize).pow(9);
    let slices = kmer_mem / max_mem + 1;
    let sz = 256 / slices + 1;

    let mut bucket_ranges = Vec::new();
    let mut start = 0;
    while start < 256 {
        bucket_ranges.push(start..start + sz);
        start += sz;
    }
    assert!(bucket_ranges[bucket_ranges.len() - 1].end >= 256);

    if bucket_ranges.len() > 1 {
        println!(
            "filter_kmers: {} sequences, {} kmers, {} passes",
            num_seqs,
            input_kmers,
            bucket_ranges.len()
        );
    }

    for bucket_range in bucket_ranges {

        let mut kmer_buckets: Vec<Vec<(K, Exts, D1)>> = Vec::new();
        for _ in 0..256 {
            kmer_buckets.push(Vec::new());
        }

        for seq in seqs {
            for (kmer, exts) in seq.iter_kmer_exts(seq_exts) {
                let (min_kmer, flip_exts) = if rc_norm {
                    let (min_kmer, flip) = kmer.min_rc_flip();
                    let flip_exts = if flip { exts.rc() } else { exts };
                    (min_kmer, flip_exts)
                } else {
                    (kmer, exts)
                };
                let bucket = filter::bucket(min_kmer);

                if bucket >= bucket_range.start && bucket < bucket_range.end {
                    kmer_buckets[bucket].push((min_kmer, flip_exts, d.clone()));
                }
            }
        }

        // info!("Validating kmers...");
        for mut kmer_vec in kmer_buckets {
            pdqsort::sort_by_key(&mut kmer_vec, |elt| elt.0);

            for (kmer, kmer_obs_iter) in &kmer_vec.into_iter().group_by(|elt| elt.0) {
                let (is_valid, exts, summary_data) = summarizer.summarize(kmer_obs_iter);
                // if report_all_kmers { all_kmers.push(kmer); }
                if is_valid {
                    valid_kmers.push(kmer);
                    valid_exts.push(exts);
                    // valid_data.push(summary_data);
                }
            }
        }
    }

    // println!("{:?}", valid_kmers.heap_size_of_children());
    // pdqsort::sort_by_key(&mut valid_kmers, |x| x.0);
    // pdqsort::sort(&mut all_kmers);
    //if report_all_kmers {
    //     filter::remove_censored_exts_sharded(stranded, &mut valid_kmers, &all_kmers);
    //}

    println!(
        "filter kmers: sequences: {}, kmers: {}, unique kmers: {}",
        num_seqs,
        input_kmers,
        // all_kmers.len(),
        valid_kmers.len()
    );

    //BoomHashMap{
    //    mph_hash: boomphf::Mphf::new(1.7, &valid_kmers, None),
    //    keys: valid_kmers,
    //    data: valid_data,
    //    exts: valid_exts,
    //}
}
