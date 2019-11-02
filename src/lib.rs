// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bincode;
extern crate bio;
extern crate boomphf;
extern crate crossbeam_utils;
extern crate debruijn;
extern crate flate2;
#[macro_use]
extern crate failure;
extern crate itertools;
extern crate rayon;
extern crate regex;
extern crate rust_htslib;
extern crate smallvec;
extern crate shardio;

#[macro_use] 
extern crate pretty_assertions;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate serde;

#[macro_use]
extern crate log;

pub mod bam;
pub mod build_index;
pub mod em;
pub mod config;
pub mod mappability;
pub mod pseudoaligner;
pub mod hla;
pub mod utils;
mod locus;

