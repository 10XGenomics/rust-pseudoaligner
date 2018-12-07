// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bincode;
extern crate bio;
extern crate boomphf;
extern crate crossbeam;
extern crate debruijn;
extern crate flate2;
extern crate failure;
extern crate itertools;
extern crate rayon;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate serde;

#[macro_use]
extern crate log;

pub mod build_index;
pub mod config;
pub mod mappability;
pub mod pseudoaligner;
pub mod utils;
