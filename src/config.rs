use debruijn::{kmer};
use std::collections::HashMap;

// main configs
pub const MEM_SIZE: usize = 1;
pub const MIN_KMERS: usize = 1;
pub const STRANDED: bool = true;
pub const REPORT_ALL_KMER: bool = false;
pub const READ_COVERAGE_THRESHOLD: usize = 32;
pub const BUCKET_SIZE_THRESHOLD: usize = 500;

pub const U8_MAX: usize = u8::max_value() as usize;
pub const U16_MAX: usize = u16::max_value() as usize;
pub const U32_MAX: usize = u32::max_value() as usize;

// Worker queue configs
pub const MAX_WORKER: usize = 6;

//DOCKS configs
pub type DocksUhs = HashMap<String, u16>;
pub const DOCKS_FILE: &'static str = "res_6_30_4_0.txt";

//pub const K: usize = Minimizer::k();
//pub const L: usize = KmerType::k();
// NOTE: Rust don't allow static from function
// Remember to change Minimizer and KmerType if
// you are changing K, L
pub type Minimizer = kmer::Kmer6;
pub const K: usize = 6;

pub type KmerType = kmer::Kmer30;
pub const L: usize = 30;
