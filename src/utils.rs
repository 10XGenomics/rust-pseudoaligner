use std;
use std::hash::Hash;
use std::fmt::Debug;
use std::marker::PhantomData;

use pdqsort;
use boomphf;
use fxhash::FxHashMap;

use itertools::Itertools;
use debruijn::{compression, filter};
use debruijn::{Dir, Kmer, Exts, Vmer, Mer};

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

// Minimal perfect hash
pub struct BoomHashMap<K: Clone + Hash + Debug, Exts, D> {
    mph_hash: boomphf::Mphf<K>,
    keys: Vec<K>,
    exts: Vec<Exts>,
    data: Vec<D>
}

/// Modified version of filter_kmer for sc_mapping
pub fn kmerize<K: Kmer, V: Vmer<K>, D: Clone + Ord>(
    seqs: Vec<(V, Exts, D)>,
) -> BoomHashMap<K, Exts, D> {

    let mut valid_kmers: FxHashMap<K, Vec<D>> = FxHashMap::default();

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs.iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
        .sum();

    println!(
        "Parse Fasta and found: {} sequences, {} input kmers",
        seqs.len(),
        input_kmers,
    );

    for &(ref seq, seq_exts, ref d) in &seqs {
        for (kmer, exts) in seq.iter_kmer_exts(seq_exts) {
            let (min_kmer, flip_exts) =  {
                let (min_kmer, flip) = kmer.min_rc_flip();
                let flip_exts = if flip { exts.rc() } else { exts };
                (min_kmer, flip_exts)
            };

            if valid_kmers.contains_key(&min_kmer){
                match valid_kmers.get_mut(&min_kmer) {
                    Some(kmer_data) => {
                        kmer_data.push(d.clone());
                        pdqsort::sort(kmer_data.as_mut_slice());
                        kmer_data.dedup();
                    },
                    None => {}
                }
            }
            else{
                valid_kmers.insert(min_kmer, vec![d.clone()]);
            }
        }
    }

    println!(
        "filter kmers: sequences: {}, Input-kmers: {}, Unique-kmers: {}",
        seqs.len(),
        input_kmers,
        valid_kmers.len()
    );

    let dummy_vec : Vec<K> = vec![K::from_ascii(b"GTTAACTTGCCGTCAGCCTTTTCTTTGACCTCTTCTTT")];
    print!("{:?}", dummy_vec);
    let phf: boomphf::Mphf<K> = boomphf::Mphf::new(1.7, &dummy_vec.clone(), None);
    BoomHashMap{
        mph_hash: phf,
        keys: Vec::new(),
        data: Vec::new(),
        exts: Vec::new(),
    }
}
