// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.


use std::fs::File;
use std::hash::Hash;
use std::fmt::Debug;
use std::boxed::Box;
use std::path::{Path};
use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter};


//use bincode;
use bincode;
use boomphf;
use serde::{Serialize};
use serde::de::DeserializeOwned;
use debruijn::graph::DebruijnGraph;
use bincode::{serialize_into, deserialize_from};

use failure::Error;
use flate2::read::MultiGzDecoder;

#[derive(Serialize, Deserialize)]
pub struct Index<K, E, D>
where K:Hash, D: Eq + Hash{
    dbg: DebruijnGraph<K, D>,
    phf: boomphf::BoomHashMap2<K, E, D>,
    eqclasses: HashMap<Vec<D>, D>,
}

impl<K:Hash, E, D> Index<K, E, D> 
where K:Hash, D: Eq + Hash{
    pub fn new(dbg: DebruijnGraph<K, D>,
               phf: boomphf::BoomHashMap2<K, E, D>,
               eqclasses: HashMap<Vec<D>, D>) -> Index<K, E, D>{
        Index{
            dbg: dbg,
            phf: phf,
            eqclasses: eqclasses,
        }
    }

    pub fn get_phf(&self) -> &boomphf::BoomHashMap2<K, E, D>{
        &self.phf
    }

    pub fn get_dbg(&self) -> &DebruijnGraph<K, D>{
        &self.dbg
    }

    pub fn get_eq_classes(&self) -> &HashMap<Vec<D>, D>{
        &self.eqclasses
    }
}

/// Open a (possibly gzipped) file into a BufReader.
pub fn open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<BufRead>, Error> {
    let r = File::open(p.as_ref())?;

    if p.as_ref().extension().unwrap() == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(32*1024, gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32*1024, r);
        Ok(Box::new(buf_reader))
    }
}



pub fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(g: &T, filename: P) -> Result<(), bincode::Error> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    serialize_into(&mut writer, &g)
}

pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(filename: P) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}
