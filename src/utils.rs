use std;
use std::io;
use std::mem;
use std::io::Write;
use std::hash::Hash;
use std::fmt::Debug;
use std::fmt::{self};
use std::marker::PhantomData;
use std::collections::VecDeque;

use pdqsort;
use boomphf;
use bit_set::BitSet;
use smallvec::SmallVec;
use itertools::Itertools;
use debruijn::graph::{BaseGraph};
use debruijn::{compression, filter};
use debruijn::{Dir, Kmer, Exts, Vmer};

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


impl<D: Ord+Debug> filter::KmerSummarizer<D, SmallVec<[D; 4]>> for CountFilterSmallInt<D> {
    fn summarize<K, F: Iterator<Item = (K, Exts, D)>>(&self, items: F) -> (bool, Exts, SmallVec<[D; 4]>) {
        let mut all_exts = Exts::empty();
        let mut out_data = SmallVec::<[D; 4]>::new();

        let mut nobs = 0;
        for (_, exts, d) in items {
            out_data.push(d);
            all_exts = all_exts.add(exts);
            nobs += 1;
        }

        out_data.sort();  out_data.dedup();

        (nobs as usize >= self.min_kmer_obs, all_exts, out_data)
    }
}


// Minimal perfect hash
pub struct BoomHashMap<K: Clone + Hash + Debug, Exts, D> {
    mphf: boomphf::Mphf<K>,
    keys: Vec<K>,
    exts: Vec<Exts>,
    data: Vec<D>
}
impl<K, Exts, D> BoomHashMap<K, Exts, D>
where K: Clone + Hash + Debug, D: Debug, Exts: Debug {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Kmer {{ id: {:?}, Exts: {:?}, Data: {:?} }} \n Lengths {}, {}, {}",
            self.keys,
            self.exts,
            self.data,
            self.keys.len(),
            self.exts.len(),
            self.data.len()
        )
    }

    fn new(mut keys_: Vec<K>, mut exts_: Vec<Exts>, mut data_: Vec<D> ) -> BoomHashMap<K, Exts, D> {
        let mphf_ = boomphf::Mphf::new(1.7, &keys_, None);
        // trick taken from :
        // https://github.com/10XDev/cellranger/blob/master/lib/rust/detect_chemistry/src/index.rs#L123
        for i in 0 .. keys_.len() {
            loop {
                let kmer_slot = mphf_.hash(&keys_[i]) as usize;
                if i == kmer_slot { break; }
                keys_.swap(i, kmer_slot);
                exts_.swap(i, kmer_slot);
                data_.swap(i, kmer_slot);
            }
        }
        BoomHashMap{
            mphf: mphf_,
            keys: keys_,
            exts: exts_,
            data: data_,
        }
    }
}

/// Read a shard and determine the valid kmers
/// Low memory implementation that should consume < 4G of temporary memory
/// To reduce memory consumption, set track_bcs to false to forget about BC lists.
#[inline(never)]
pub fn filter_kmers_with_mphf<K: Kmer, V: Vmer<K>, D1: Clone, DS, S: filter::KmerSummarizer<D1, DS>>(
    seqs: Vec<(V, Exts, D1)>,
    summarizer: S,
    stranded: bool,
    //report_all_kmers: bool,
    memory_size: usize,
) -> BoomHashMap<K, Exts, DS>
where DS: std::fmt::Debug{

    let rc_norm = !stranded;

    // let mut all_kmers = Vec::new();
    let mut valid_kmers = Vec::new();
    let mut valid_exts = Vec::new();
    let mut valid_data = Vec::new();

    let num_seqs: usize = seqs.len();

    // Estimate memory consumed by Kmer vectors, and set iteration count appropriately
    let input_kmers: usize = seqs.iter()
        .map(|&(ref vmer, _, _)| vmer.len().saturating_sub(K::k() - 1))
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
            "\nfilter_kmers: {} sequences, {} kmers, {} passes",
            num_seqs,
            input_kmers,
            bucket_ranges.len()
        );
    }

    let mut pass_counter = 0;
    for bucket_range in bucket_ranges {

        let mut kmer_buckets: Vec<Vec<(K, Exts, D1)>> = Vec::new();
        pass_counter += 1;
        print!("\r Performing {} pass", pass_counter);
        io::stdout().flush().ok().expect("Could not flush stdout");

        for _ in 0..256 {
            kmer_buckets.push(Vec::new());
        }

        for &(ref seq, seq_exts, ref d) in &seqs {
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
                    valid_data.push(summary_data);
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
        "\nfilter kmers: sequences: {}, kmers: {}, unique kmers: {}",
        num_seqs,
        input_kmers,
        // all_kmers.len(),
        valid_kmers.len()
    );

    BoomHashMap::new(valid_kmers, valid_exts, valid_data)
}

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

//////////////////////////////
// Compress from Hash a new Struct
//////////////////////////////
/// Generate a compressed DeBruijn graph from hash_index
struct CompressFromHash<'a, K: 'a + Kmer, D: 'a, S: compression::CompressionSpec<D>> {
    stranded: bool,
    k: PhantomData<K>,
    d: PhantomData<D>,
    spec: S,
    available_kmers: BitSet,
    index: &'a BoomHashMap<K, Exts, D>,
}

/// Compression of paths in Debruijn graph
impl<'a, K: Kmer, D: Clone + Debug, S: compression::CompressionSpec<D>> CompressFromHash<'a, K, D, S> {
    fn get_kmer_data(&'a self, kmer: &K) -> (Exts, &'a D) {
        let pos = match self.get_kmer_id(kmer).ok() {
            Some(i) => i,
            None => panic!("couldn't find kmer {:?}", kmer),
        };

        (self.index.exts[pos], &self.index.data[pos])
    }

    fn get_kmer_id(&self, kmer: &K) -> Result<usize, usize> {
        self.index.mphf.try_hash(kmer).map_or(None, |v| Some(v as usize)).ok_or_else(|| 0 as usize)
    }

    /// Attempt to extend kmer v in direction dir. Return:
    ///  - Unique(nextKmer, nextDir) if a single unique extension
    ///    is possible.  nextDir indicates the direction to extend nextMker
    ///    to preserve the direction of the extension.
    /// - Term(ext) no unique extension possible, indicating the extensions at this end of the line
    fn try_extend_kmer(&self, kmer: K, dir: Dir) -> compression::ExtMode<K> {

        // metadata of start kmer
        let (exts, ref kmer_data) = self.get_kmer_data(&kmer);

        if exts.num_ext_dir(dir) != 1 || (!self.stranded && kmer.is_palindrome()) {
            compression::ExtMode::Terminal(exts.single_dir(dir))
        } else {
            // Get the next kmer
            let ext_base = exts.get_unique_extension(dir).expect("should be unique");

            let mut next_kmer = kmer.extend(ext_base, dir);

            let mut do_flip = false;

            if !self.stranded {
                let flip_rc = next_kmer.min_rc_flip();
                do_flip = flip_rc.1;
                next_kmer = flip_rc.0;
            }

            let next_dir = dir.cond_flip(do_flip);
            let is_palindrome = !self.stranded && next_kmer.is_palindrome();


            // We can include this kmer in the line if:
            // a) it exists in the partition and is still unused
            // b) the kmer we go to has a unique extension back in our direction

            // Check condition a)
            match self.get_kmer_id(&next_kmer) {
                Ok(id) if self.available_kmers.contains(id) => (),

                // This kmer isn't in this partition, or we've already used it
                _ => return compression::ExtMode::Terminal(exts.single_dir(dir))
            }

            // Check condition b)
            // Direction we're approaching the new kmer from
            let new_incoming_dir = dir.flip().cond_flip(do_flip);
            let next_kmer_r = self.get_kmer_data(&next_kmer);
            let (next_kmer_exts, ref next_kmer_data) = next_kmer_r;
            let incoming_count = next_kmer_exts.num_ext_dir(new_incoming_dir);
            let outgoing_exts = next_kmer_exts.single_dir(new_incoming_dir.flip());

            // Test if the spec let's us combine these into the same path
            let can_join = self.spec.join_test(kmer_data, next_kmer_data);

            if incoming_count == 0 && !is_palindrome {
                panic!("unreachable");
            } else if can_join && incoming_count == 1 && !is_palindrome {
                // We have a unique path to next_kmer -- include it
                compression::ExtMode::Unique(next_kmer, next_dir, outgoing_exts)
            } else {
                // there's more than one path
                // into the target kmer - don't include it
                compression::ExtMode::Terminal(exts.single_dir(dir))
            }
        }
    }


    /// Build the maximal line starting at kmer in direction dir, at most max_dist long.
    /// Also return the extensions at the end of this line.
    /// Sub-lines break if their extensions are not available in this shard
    #[inline(never)]
    fn extend_kmer(&mut self, kmer: K, start_dir: Dir, path: &mut Vec<(K, Dir)>) -> Exts {

        let mut current_dir = start_dir;
        let mut current_kmer = kmer;
        path.clear();

        let final_exts: Exts; // must get set below

        let id = self.get_kmer_id(&kmer).expect("should have this kmer");
        let _ = self.available_kmers.remove(id);

        loop {
            let ext_result = self.try_extend_kmer(current_kmer, current_dir);

            match ext_result {
                compression::ExtMode::Unique(next_kmer, next_dir, _) => {
                    path.push((next_kmer, next_dir));
                    let next_id = self.get_kmer_id(&next_kmer).expect("should have this kmer");
                    self.available_kmers.remove(next_id);
                    current_kmer = next_kmer;
                    current_dir = next_dir;
                }
                compression::ExtMode::Terminal(ext) => {
                    final_exts = ext;
                    break;
                }
            }
        }

        final_exts
    }


    /// Build the edge surrounding a kmer
    #[inline(never)]
    fn build_node(
        &mut self,
        seed: K,
        path: &mut Vec<(K, Dir)>,
        edge_seq: &mut VecDeque<u8>,
    ) -> (Exts, D) {

        edge_seq.clear();
        for i in 0..K::k() {
            edge_seq.push_back(seed.get(i));
        }

        let mut node_data = self.get_kmer_data(&seed).1.clone();

        let l_ext = self.extend_kmer(seed, Dir::Left, path);

        // Add on the left path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer,
                Dir::Right => next_kmer.rc(),
            };

            edge_seq.push_front(kmer.get(0));

            // Reduce the data object
            let (_, ref kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let left_extend = match path.last() {
            None => l_ext,
            Some(&(_, Dir::Left)) => l_ext,
            Some(&(_, Dir::Right)) => l_ext.complement(),
        };

        let r_ext = self.extend_kmer(seed, Dir::Right, path);

        // Add on the right path
        for &(next_kmer, dir) in path.iter() {
            let kmer = match dir {
                Dir::Left => next_kmer.rc(),
                Dir::Right => next_kmer,
            };

            edge_seq.push_back(kmer.get(K::k() - 1));

            let (_, ref kmer_data) = self.get_kmer_data(&next_kmer);
            node_data = self.spec.reduce(node_data, kmer_data)
        }

        let right_extend = match path.last() {
            None => r_ext,
            Some(&(_, Dir::Left)) => r_ext.complement(),
            Some(&(_, Dir::Right)) => r_ext,
        };

        (Exts::from_single_dirs(left_extend, right_extend), node_data)
    }

    /// Compress a set of kmers and their extensions and metadata into a base DeBruijn graph.
    #[inline(never)]
    pub fn compress_kmers(
        stranded: bool,
        spec: S,
        index: &BoomHashMap<K, Exts, D>,
    ) -> BaseGraph<K, D> {

        let n_kmers = index.keys.len();
        let mut available_kmers = BitSet::with_capacity(n_kmers);
        for i in 0..n_kmers {
            available_kmers.insert(i);
        }

        let mut comp = CompressFromHash {
            stranded: stranded,
            spec: spec,
            k: PhantomData,
            d: PhantomData,
            available_kmers: available_kmers,
            index: index,
        };

        // Path-compressed De Bruijn graph will be created here
        let mut graph = BaseGraph::new(stranded);

        // Paths will be get assembled here
        let mut path_buf = Vec::new();

        // Node sequences will get assembled here
        let mut edge_seq_buf = VecDeque::new();

        for kmer_counter in 0..n_kmers {
            let start_kmer = index.keys[kmer_counter].clone();
            if comp.available_kmers.contains(kmer_counter) {
                let (node_exts, node_data) =
                    comp.build_node(start_kmer, &mut path_buf, &mut edge_seq_buf);
                graph.add(&edge_seq_buf, node_exts, node_data);
            }
        }

        graph
    }
}

/// Take a BoomHash Object and build a compressed DeBruijn graph.
#[inline(never)]
pub fn compress_kmers_with_hash<K: Kmer, D: Clone + Debug, S: compression::CompressionSpec<D>>(
    stranded: bool,
    spec: S,
    index: &BoomHashMap<K, Exts, D>,
) -> BaseGraph<K, D> {
    CompressFromHash::<K, D, S>::compress_kmers(stranded, spec, index)
}
