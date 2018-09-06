// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bincode;
extern crate bio;
extern crate boomphf;
extern crate crossbeam;
extern crate csv;
extern crate debruijn;
extern crate docopt;
extern crate failure;
extern crate flate2;
extern crate itertools;
extern crate num;
extern crate pretty_env_logger;
extern crate rayon;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate log;

#[macro_use]
extern crate serde;

pub mod build_index;
pub mod config;
pub mod pseudoaligner;
pub mod utils;
