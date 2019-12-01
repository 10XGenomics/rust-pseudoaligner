// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use lazy_static::lazy_static;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::config::{KmerType, MEM_SIZE, REPORT_ALL_KMER, STRANDED};
use boomphf::hashmap::{BoomHashMap2, NoKeyBoomHashMap};
use debruijn;
use debruijn::compression::*;
use debruijn::dna_string::{DnaString, DnaStringSlice};
use debruijn::filter::*;
use debruijn::graph::*;
use debruijn::*;

use crate::config::{MIN_KMERS, U32_MAX};
use crate::equiv_classes::{CountFilterEqClass, EqClassIdType};
use crate::pseudoaligner::Pseudoaligner;
use boomphf;
use failure::Error;
use log::info;
use rayon::prelude::*;
use rayon::{self, ThreadPool};

const MIN_SHARD_SEQUENCES: usize = 2000;

pub fn build_index<K: Kmer + Sync + Send>(
    seqs: &[DnaString],
    tx_names: &Vec<String>,
    tx_gene_map: &HashMap<String, String>,
    num_threads: usize,
) -> Result<Pseudoaligner<K>, Error> {
    // Thread pool Configuration for calling BOOMphf
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build()?;

    if seqs.len() >= U32_MAX {
        panic!("Too many ({}) sequences to handle.", seqs.len());
    }

    info!("Sharding sequences...");

    let mut buckets: Vec<_> = seqs
        .iter()
        .enumerate()
        .flat_map(|(id, seq)| partition_contigs::<KmerType>(seq, id as u32))
        .collect();

    pool.install(|| {
        buckets.par_sort_unstable_by_key(|x| x.0);
    });
    info!("Got {} sequence chunks", buckets.len());

    let summarizer = Arc::new(CountFilterEqClass::new(MIN_KMERS));
    let sequence_shards = group_by_slices(&buckets, |x| x.0, MIN_SHARD_SEQUENCES);

    info!("Assembling {} shards...", sequence_shards.len());

    let shard_dbgs = pool.install(|| {
        let mut shard_dbgs = Vec::with_capacity(sequence_shards.len());
        sequence_shards
            .into_par_iter()
            .into_par_iter()
            .map_with(summarizer.clone(), |s, strings| {
                assemble_shard::<K>(strings, s)
            })
            .collect_into_vec(&mut shard_dbgs);

        shard_dbgs
    });

    info!("Done dBG construction of shards");
    info!("Starting merging disjoint graphs");

    let dbg = merge_shard_dbgs(shard_dbgs);
    info!("Graph merge complete");

    let eq_classes = summarizer.get_eq_classes();

    info!("Indexing de Bruijn graph");
    let dbg_index = make_dbg_index(&dbg, &pool, num_threads);

    Ok(Pseudoaligner::new(
        dbg,
        eq_classes,
        dbg_index,
        tx_names.clone(),
        tx_gene_map.clone(),
    ))
}

// Manually compute the equivalence class of each kmer, and make sure
// it matches that equivalence class for that kmer inside the DBG.
#[inline(never)]
pub fn validate_dbg<K: Kmer + Sync + Send>(seqs: &[DnaString], al: &Pseudoaligner<K>) {
    let mut eqclasses = HashMap::<K, Vec<u32>>::new();

    // compute the equivalence class of each kmer
    for (i, s) in seqs.iter().enumerate() {
        for k in s.iter_kmers::<K>() {
            let eq = eqclasses.entry(k).or_insert(vec![]);
            eq.push(i as u32)
        }
    }

    // check that the equivalence class of the kmer inside the graph matches the naive version
    for (k, mut test_eqclass) in eqclasses {
        test_eqclass.dedup();

        if test_eqclass.len() > 5000 {
            println!("kmer: {:?}, eqclass.len(): {}", k, test_eqclass.len());
        }

        let (node_id, _) = al.dbg_index.get(&k).unwrap();

        let eq_class = al.dbg.get_node(*node_id as usize).data();
        let dbg_eqclass = &al.eq_classes[*eq_class as usize];

        let mut dbg_eq_clone = dbg_eqclass.clone();
        dbg_eq_clone.dedup();

        if &dbg_eq_clone != dbg_eqclass {
            println!(
                "dbg eq class not unique: eqclass_id: {}, node: {}",
                eq_class, node_id
            );
            assert_eq!(&dbg_eq_clone, dbg_eqclass);
        }

        assert_eq!(&test_eqclass, dbg_eqclass);
    }

    // check that each read sequence aligns cleanly
    for (i, s) in seqs.iter().enumerate() {
        let i = i as u32;

        // transcripts shorter than k can't be mapped
        if s.len() < K::k() {
            continue;
        }

        let (eqclass, bases_aligned) = al.map_read(s).unwrap();
        assert_eq!(s.len(), bases_aligned);

        if eqclass.len() > 1 {
            assert!(eqclass.contains(&i));

            // identical strings
            if eqclass.len() == 2 && seqs[eqclass[0] as usize] == seqs[eqclass[1] as usize] {
                continue;
            }

            // if the sequences aren't identical, the current string must be shortest, or
            // the set of nodes visited by the input string must be a subset of the other sequences
            // in the equivalence class.
            let shortest = eqclass
                .iter()
                .map(|x| seqs[*x as usize].len())
                .min()
                .unwrap();

            if s.len() != shortest {
                let mut path_buf: Vec<usize> = Vec::new();

                use std::iter::FromIterator;
                al.map_read_to_nodes(s, &mut path_buf).unwrap();
                let my_nodes: HashSet<usize> = HashSet::from_iter(path_buf.iter().cloned());

                println!("eqclass: {:?}", eqclass);
                for i in &eqclass {
                    println!(
                        "{}: {}, len:{}",
                        i,
                        al.tx_names[*i as usize],
                        seqs[*i as usize].len()
                    );
                    println!("{:?}", seqs[*i as usize]);

                    let r = al
                        .map_read_to_nodes(&seqs[*i as usize], &mut path_buf)
                        .unwrap();
                    let other_nodes = HashSet::from_iter(path_buf.iter().cloned());

                    println!("r: {:?}", r);
                    println!("{:?}", path_buf);

                    assert!(my_nodes.is_subset(&other_nodes));
                    println!("---");
                }
            }

        // debugging
        // println!("--- dup on {}", i);
        // for e in &eqclass {
        //     println!("{}: {}", e, al.tx_names[*e as usize]);
        //     println!("{:?}", seqs[*e as usize]);
        // }
        } else {
            assert_eq!(eqclass, vec![i]);
        }
    }
}

type PmerType = debruijn::kmer::Kmer6;

lazy_static! {
    static ref PERM: Vec<usize> = {
        let maxp = 1 << (2 * PmerType::k());
        let mut sorted_kmers: Vec<PmerType> = Vec::with_capacity(maxp);
        for i in 0..maxp {
            let kmer = Kmer::from_u64(i as u64);
            sorted_kmers.push(kmer);
        }
        sorted_kmers.sort_by_key(count_a_t_bases);

        let mut permutation = Vec::new();
        permutation.resize(maxp, 0);

        for (sort_pos, kmer) in sorted_kmers.into_iter().enumerate() {
            permutation[kmer.to_u64() as usize] = sort_pos;
        }

        permutation
    };
}

/// Count the number of A/T bases in a kmer
fn count_a_t_bases<K: Kmer>(kmer: &K) -> usize {
    let mut count = 0;
    for i in 0..K::k() {
        if kmer.get(i) == b'A' || kmer.get(i) == b'T' {
            count += 1;
        }
    }
    count
}

fn partition_contigs<'a, K: Kmer>(
    contig: &'a DnaString,
    contig_id: u32,
) -> Vec<(u16, u32, DnaStringSlice<'a>, Exts)> {
    // One FASTA entry possibly broken into multiple contigs
    // based on the location of `N` int he sequence.

    let mut bucket_slices = Vec::new();

    if contig.len() >= K::k() {
        let msps = debruijn::msp::simple_scan::<_, PmerType>(K::k(), contig, &PERM, false);
        for msp in msps {
            let bucket_id = msp.bucket();
            let slice = contig.slice(msp.start(), msp.end());
            let exts = Exts::from_dna_string(contig, msp.start(), msp.len());
            bucket_slices.push((bucket_id, contig_id, slice, exts));
        }
    }

    bucket_slices
}

fn assemble_shard<K: Kmer>(
    shard_data: &[(u16, u32, DnaStringSlice, Exts)],
    summarizer: &Arc<CountFilterEqClass<u32>>,
) -> BaseGraph<K, EqClassIdType> {
    let filter_input: Vec<_> = shard_data
        .into_iter()
        .cloned()
        .map(|(_, seqid, string, exts)| (string, exts, seqid))
        .collect();

    let (phf, _): (BoomHashMap2<K, Exts, EqClassIdType>, _) = filter_kmers(
        &filter_input,
        summarizer,
        STRANDED,
        REPORT_ALL_KMER,
        MEM_SIZE,
    );

    compress_kmers_with_hash(STRANDED, &ScmapCompress::new(), &phf)
}

fn merge_shard_dbgs<K: Kmer + Sync + Send>(
    uncompressed_dbgs: Vec<BaseGraph<K, EqClassIdType>>,
) -> DebruijnGraph<K, EqClassIdType> {
    let combined_graph = BaseGraph::combine(uncompressed_dbgs.into_iter()).finish();
    compress_graph(STRANDED, &ScmapCompress::new(), combined_graph, None)
}

#[inline(never)]
fn make_dbg_index<K: Kmer + Sync + Send>(
    dbg: &DebruijnGraph<K, EqClassIdType>,
    pool: &ThreadPool,
    num_threads: usize,
) -> NoKeyBoomHashMap<K, (u32, u32)> {
    let mut total_kmers = 0;
    let kmer_length = K::k();
    for node in dbg.iter_nodes() {
        total_kmers += node.len() - kmer_length + 1;
    }

    println!("Total {:?} kmers to process in dbg", total_kmers);
    println!("Making mphf of kmers");
    let mphf =
        //boomphf::Mphf::from_chunked_iterator(1.7, dbg, total_kmers);
        boomphf::Mphf::from_chunked_iterator_parallel(1.7, dbg, None, total_kmers, num_threads);

    println!("Assigning offsets to kmers");
    let mut node_and_offsets = Vec::with_capacity(total_kmers);
    node_and_offsets.resize(total_kmers, (U32_MAX as u32, U32_MAX as u32));

    let scatter = crate::scatter::ScatterToVec::new(&mut node_and_offsets[..]);

    pool.install(|| {
        (0..dbg.len()).into_par_iter().for_each_init(
            // Each thread gets a scatter handle to write values from
            || scatter.handle(),
            |handle, node_id| {
                let node = dbg.get_node_kmer(node_id);

                for (offset, kmer) in node.into_iter().enumerate() {
                    let index = mphf.try_hash(&kmer).expect("can't find kmer in DBG graph!");
                    handle.write(index as usize, (node_id as u32, offset as u32));
                }
            },
        );
    });

    boomphf::hashmap::NoKeyBoomHashMap::new_with_mphf(mphf, node_and_offsets)
}

/// Split the slice `data` into subslices of size at least
/// `min_size`, while ensuring that consecutive runs of
/// items with the same key as defined by the key function `f` are
/// in the same subslice.
fn group_by_slices<T, K: PartialEq, F: Fn(&T) -> K>(
    data: &[T],
    f: F,
    min_size: usize,
) -> Vec<&[T]> {
    let mut slice_start = 0;
    let mut result = Vec::new();
    for i in 1..data.len() {
        if !(f(&data[i - 1]) == f(&data[i])) && (i - slice_start) > min_size {
            result.push(&data[slice_start..i]);
            slice_start = i;
        }
    }
    if slice_start < data.len() {
        result.push(&data[slice_start..]);
    }
    result
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::config;
    use crate::utils;
    use bio::io::fasta;
    use failure::ResultExt;
    use proptest::collection::vec;
    use proptest::prelude::*;
    use proptest::proptest;

    proptest! {
        #![proptest_config(ProptestConfig { cases: 2000, .. ProptestConfig::default()})]
        #[test]
        fn group_by_slices_test(
            //v: Vec<u16>,
            v in vec(0..100usize, 0..5000usize),
            min_sz in 1..200usize
        ) {
            let res =  group_by_slices(&v, |v| v.clone(), min_sz);
            let total_len: usize = res.iter().map(|x| x.len()).sum();
            prop_assert_eq!(v.len(), total_len);

            for i in 1 .. res.len() {
                prop_assert!(res[i-1].len() >= min_sz);
            }

            for i in 1 .. res.len() {
                let prev = res[i-1];
                let next = res[i];
                prop_assert!(prev[prev.len() - 1] != next[0]);
            }
        }
    }

    #[test]
    fn test_gencode_small_build() -> Result<(), Error> {
        let fasta = fasta::Reader::from_file("test/gencode_small.fa")?;
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(fasta)?;
        let index = build_index::<config::KmerType>(&seqs, &tx_names, &tx_gene_map, 2)?;
        validate_dbg(&seqs, &index);
        Ok(())
    }

    #[cfg_attr(feature = "slow_tests", test)]
    fn test_gencode_full_build() -> Result<(), Error> {
        let msg = "For full txome indexing test, download from ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.transcripts.fa.gz, un-gzip and place in test/gencode.v28.transcripts.fa";
        let fasta = fasta::Reader::from_file("test/gencode.v28.transcripts.fa").context(msg)?;
        let (seqs, tx_names, tx_gene_map) = utils::read_transcripts(fasta)?;
        let index = build_index::<config::KmerType>(&seqs, &tx_names, &tx_gene_map, 2)?;
        validate_dbg(&seqs, &index);
        Ok(())
    }
}
