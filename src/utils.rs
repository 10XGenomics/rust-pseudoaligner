// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.

use std::fs;
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
use debruijn::{Kmer, Vmer};
use serde::{Serialize};
use serde::de::DeserializeOwned;
use debruijn::graph::DebruijnGraph;
use debruijn::filter::EqClassIdType;
use bincode::{serialize_into, deserialize_from};

use num;
use num::NumCast;
use failure::Error;
use flate2::read::MultiGzDecoder;

#[derive(Serialize, Deserialize)]
pub struct Index<K, D>
where K:Hash + Serialize, D: Eq + Hash + Serialize {
    eqclasses: Vec<Vec<D>>,
    dbg: DebruijnGraph<K, EqClassIdType>,
    left_phf: boomphf::BoomHashMap<K, EqClassIdType>,
    right_phf: boomphf::BoomHashMap<K, EqClassIdType>,
}

impl<K, D> Index<K, D>
where K:Hash + Serialize + Kmer + Send + Sync + DeserializeOwned,
      D: Clone + Debug + Eq + Hash + Serialize + NumCast + DeserializeOwned{
    pub fn dump(dbg: DebruijnGraph<K, EqClassIdType>,
                eqclasses: HashMap<Vec<D>, EqClassIdType>,
                index_path: &str) {

        let hash_len = eqclasses.len();
        let mut eqclasses_vec: Vec<Vec<D>> = Vec::new();
        eqclasses_vec.resize(hash_len, Vec::<D>::new());

        // fill in eqclasses into vector
        for (eqclass, index) in eqclasses.into_iter() {
            eqclasses_vec[index as usize] = eqclass;
        }
        info!("Found {} Equivalence classes", eqclasses_vec.len());

        info!("Dumping index into folder: {:?}", index_path);
        match fs::create_dir(index_path) {
            Err(err) => warn!("{:?}", err),
            Ok(()) => info!("Creating folder {:?}", index_path),
        }

        let eqclass_file_name = index_path.to_owned() + "/eq_classes.bin";
        write_obj(&eqclasses_vec, eqclass_file_name).expect("Can't dump classes");

        let dbg_file_name = index_path.to_owned() + "/dbg.bin";
        write_obj(&dbg, dbg_file_name).expect("Can't dump debruijn graph");

        let mut left_kmers = Vec::new();
        //let mut right_kmers = Vec::new();
        let mut node_ids = Vec::new();

        for node in dbg.iter_nodes() {
            //let (left_kmer, right_kmer): (K, K) = node.sequence().both_term_kmer();
            let left_kmer: K = node.sequence().first_kmer();
            left_kmers.push(left_kmer);
            //right_kmers.push(right_kmer);
            node_ids.push(node.node_id);
        }

        { // Clearing after writing
            let left_phf = boomphf::BoomHashMap::new(left_kmers, node_ids);
            let left_phf_file_name = index_path.to_owned() + "/left_phf.bin";
            write_obj(&left_phf, left_phf_file_name).expect("Can't dump left phf");
        }

        //let right_phf = boomphf::BoomHashMap::new(right_kmers, node_ids);
        //let right_phf_file_name = index_path.to_owned() + "/right_phf.bin";
        //write_obj(&right_phf, right_phf_file_name).expect("Can't dump right phf");

        let data_type: D = num::cast(0).unwrap();
        write_obj(&data_type, index_path.to_owned() + "/type.bin").expect("Can't dump data type");
    }

    pub fn read(index_path: &str) -> Index<K, D> {

        match fs::read_dir(index_path) {
            Err(_) => panic!("{:?} directory not found", index_path),
            Ok(_) => info!("Reading index from folder: {:?}", index_path),
        }
        let _data_type: D;
        _data_type = read_obj(index_path.to_owned() + "/type.bin").expect("Can't read data type");

        let eqclass_file_name = index_path.to_owned() + "/eq_classes.bin";
        let eq_classes: Vec<Vec<D>> = read_obj(eqclass_file_name).expect("Can't read classes");

        let dbg_file_name = index_path.to_owned() + "/dbg.bin";
        let dbg = read_obj(dbg_file_name).expect("Can't read debruijn graph");

        let left_phf_file_name = index_path.to_owned() + "/left_phf.bin";
        let left_phf = read_obj(left_phf_file_name).expect("Can't read left phf");

        let right_phf_file_name = index_path.to_owned() + "/right_phf.bin";
        let right_phf = read_obj(right_phf_file_name).expect("Can't read right phf");

        Index{
            eqclasses: eq_classes,
            dbg: dbg,
            left_phf: left_phf,
            right_phf: right_phf,
        }
    }

    pub fn _get_left_phf(&self) -> &boomphf::BoomHashMap<K, EqClassIdType>{
        &self.left_phf
    }

    pub fn _get_right_phf(&self) -> &boomphf::BoomHashMap<K, EqClassIdType>{
        &self.right_phf
    }

    pub fn get_dbg(&self) -> &DebruijnGraph<K, EqClassIdType>{
        &self.dbg
    }

    pub fn get_eq_classes(&self) -> &Vec<Vec<D>>{
        &self.eqclasses
    }
}

/// Open a (possibly gzipped) file into a BufReader.
fn open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<BufRead>, Error> {
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



fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(g: &T, filename: P) -> Result<(), bincode::Error> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    serialize_into(&mut writer, &g)
}

fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(filename: P) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}
