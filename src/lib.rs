// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

extern crate bincode;
extern crate bio;
extern crate boomphf;
extern crate debruijn;
extern crate failure;
extern crate rayon;

#[macro_use]
extern crate lazy_static;

#[macro_use]
extern crate serde;

pub mod build_index;
pub mod config;
pub mod pseudoaligner;
pub mod utils;
