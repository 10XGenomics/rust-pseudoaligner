// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use debruijn::kmer;

// main configs
pub const MEM_SIZE: usize = 1;
pub const MIN_KMERS: usize = 1;
pub const STRANDED: bool = true;
pub const REPORT_ALL_KMER: bool = false;
pub const READ_COVERAGE_THRESHOLD: usize = 32;

pub const U32_MAX: usize = u32::max_value() as usize;

// Worker queue configs
pub const MAX_WORKER: usize = 2;

pub type KmerType = kmer::Kmer24;
