// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.
use std::fs;
use std::mem;
use std::fs::File;
use std::io::Write;
use std::hash::Hash;
use std::fmt::Debug;
use std::boxed::Box;
use std::path::{Path};
use std::sync::{Arc, Mutex};
use std::sync::mpsc::sync_channel;
use std::io::{BufRead, BufReader, BufWriter};

//use bincode;
use bincode;
use crossbeam;
use debruijn::Kmer;
use serde::{Serialize};
use serde::de::DeserializeOwned;
use debruijn::graph::DebruijnGraph;
use debruijn::filter::EqClassIdType;
use bincode::{serialize_into, deserialize_from};

use boomphf;
use boomphf::hashmap::{NoKeyBoomHashMap2};

use failure::Error;
use flate2::read::MultiGzDecoder;

use config::MAX_WORKER;

#[derive(Serialize, Deserialize, Debug)]
pub struct Index<K, D>
where K:Hash + Serialize, D: Eq + Hash + Serialize {
    eqclasses: Vec<Vec<D>>,
    dbg: DebruijnGraph<K, EqClassIdType>,
    phf: NoKeyBoomHashMap2<K, usize, u32>,
}

impl<K, D> Index<K, D>
where K:Hash + Serialize + Kmer + Send + Sync + DeserializeOwned + Send + Sync,
      D: Clone + Debug + Eq + Hash + Serialize + DeserializeOwned{
    pub fn dump(dbg: DebruijnGraph<K, EqClassIdType>,
                gene_order: Vec<String>,
                eqclasses: Vec<Vec<D>>,
                index_path: &str) {

        info!("Dumping index into folder: {:?}", index_path);
        match fs::create_dir(index_path) {
            Err(err) => warn!("{:?}", err),
            Ok(()) => info!("Creating folder {:?}", index_path),
        }

        let data_type: usize = mem::size_of::<D>();
        write_obj(&data_type, index_path.to_owned() + "/type.bin").expect("Can't dump data type");

        let eqclass_file_name = index_path.to_owned() + "/eq_classes.bin";
        write_obj(&eqclasses, eqclass_file_name).expect("Can't dump classes");

        info!("Found {} Equivalence classes", eqclasses.len());

        let genes_file_name = index_path.to_owned() + "/genes.txt";
        let mut file_handle = File::create(genes_file_name).expect("Unable to create file");
        for i in gene_order{
            write!(file_handle, "{}\n", i).expect("can't write gene names");
        }

        let dbg_file_name = index_path.to_owned() + "/dbg.bin";
        write_obj(&dbg, dbg_file_name).expect("Can't dump debruijn graph");

        let mut total_kmers = 0;
        let kmer_length = K::k();
        for node in dbg.iter_nodes() {
            total_kmers += node.len()-kmer_length+1;
        }

        info!("Total {:?} kmers to process in dbg", total_kmers);
        let mphf = boomphf::Mphf::new_parallel_with_keys(1.7, &dbg, None,
                                                         total_kmers,
                                                         MAX_WORKER);

        let phf_file_name = index_path.to_owned() + "/phf.bin";
        write_obj(&mphf, phf_file_name).expect("Can't dump phf");

        const IO_THREADS: usize = 3;
        info!("Done Creating Minimum Perfect Hash");
        info!("Aligning positions with MPHF");
        {
            let mut node_ids = vec![0; total_kmers];
            let mphf_ref = &mphf;
            let (tx, rx) = sync_channel(IO_THREADS);

            let work_queue = Arc::new(Mutex::new(dbg.into_iter()));

            crossbeam::scope(|scope| {

                for _ in 0 .. IO_THREADS {
                    let work_queue = work_queue.clone();
                    let tx = tx.clone();

                    scope.spawn(move || {
                        loop {

                            let node =
                                match work_queue.lock().unwrap().next() {
                                    None => break,
                                    Some(val) => val,
                                };

                            let node_id = node.node_id;
                            let mut indices: Vec<u64> = Vec::new();

                            for kmer in node {
                                let index = match mphf_ref.try_hash(&kmer) {
                                    None => panic!("can't find in hash"),
                                    Some(index) => index.clone(),
                                };
                                indices.push(index);
                            }

                            tx.send((indices, node_id)).unwrap();

                        } //end-loop
                    }); //end-scope
                } //end-threads-for

                drop(tx);
                for (data, node_id) in rx.iter(){
                    for index in data {
                        node_ids[index as usize] = node_id;
                    }
                }
            }); //end-crossbeam

            info!("Starting Dumping positions");
            let pos_file_name = index_path.to_owned() + "/positions.bin";
            write_obj(&node_ids, pos_file_name).expect("Can't dump positions");
        } //end-positions scope

        info!("Aligning offsets with MPHF");
        {
            let mut offsets: Vec<u32> = vec![0; total_kmers];
            let mphf_ref = &mphf;
            let (tx, rx) = sync_channel(IO_THREADS);

            let work_queue = Arc::new(Mutex::new(dbg.into_iter()));

            crossbeam::scope(|scope| {

                for _ in 0 .. IO_THREADS {
                    let work_queue = work_queue.clone();
                    let tx = tx.clone();

                    scope.spawn(move || {
                        loop {

                            let node =
                                match work_queue.lock().unwrap().next() {
                                    None => break,
                                    Some(val) => val,
                                };

                            let mut indices: Vec<u64> = Vec::new();

                            for kmer in node {
                                let index = match mphf_ref.try_hash(&kmer) {
                                    None => panic!("can't find in hash"),
                                    Some(index) => index.clone(),
                                };
                                indices.push(index);
                            }

                            tx.send(indices).unwrap();

                        } //end-loop
                    }); //end-scope
                } //end-threads-for

                drop(tx);
                for data in rx.iter(){
                    for (offset, index) in data.iter().enumerate() {
                        offsets[*index as usize] = offset as u32;
                    }
                }
            }); //end-crossbeam

            info!("Starting Dumping offsets");
            let off_file_name = index_path.to_owned() + "/offsets.bin";
            write_obj(&offsets, off_file_name).expect("Can't dump offsets");
        } //end -offset scope
    }

    pub fn read(index_path: &str) -> Index<K, D> {

        match fs::read_dir(index_path) {
            Err(_) => panic!("{:?} directory not found", index_path),
            Ok(_) => info!("Reading index from folder: {:?}", index_path),
        }

        let eqclass_file_name = index_path.to_owned() + "/eq_classes.bin";
        let eq_classes: Vec<Vec<D>> = read_obj(eqclass_file_name)
            .expect("Can't read classes");

        let dbg_file_name = index_path.to_owned() + "/dbg.bin";
        let dbg = read_obj(dbg_file_name).expect("Can't read debruijn graph");

        let phf_file_name = index_path.to_owned() + "/phf.bin";
        let phf = read_obj(phf_file_name).expect("Can't read phf");

        let pos_file_name = index_path.to_owned() + "/positions.bin";
        let positions = read_obj(pos_file_name).expect("Can't read pos");

        let off_file_name = index_path.to_owned() + "/offsets.bin";
        let offsets = read_obj(off_file_name).expect("Can't read offsets");

        let phf = NoKeyBoomHashMap2::new_with_mphf( phf, positions,
                                                             offsets );
        Index{
            eqclasses: eq_classes,
            dbg: dbg,
            phf: phf,
        }
    }

    pub fn get_phf(&self) -> &NoKeyBoomHashMap2<K, usize, u32>{
        &self.phf
    }

    pub fn get_dbg(&self) -> &DebruijnGraph<K, EqClassIdType>{
        &self.dbg
    }

    pub fn get_eq_classes(&self) -> &Vec<Vec<D>>{
        &self.eqclasses
    }
}

pub fn get_data_type(index_path: &str) -> usize {
    match fs::read_dir(index_path) {
        Err(_) => panic!("{:?} directory not found", index_path),
        Ok(_) => info!("Reading index from folder: {:?}", index_path),
    }

    let data_type_file_name = index_path.to_owned() + "/type.bin";
    let data_type: usize = read_obj(data_type_file_name).expect("Can't read data type");
    data_type
}

/// Open a (possibly gzipped) file into a BufReader.
fn _open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<BufRead>, Error> {
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

pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(filename: P) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}
