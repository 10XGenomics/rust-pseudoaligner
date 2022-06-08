// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utility methods.
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter};
use std::path::Path;
use std::sync::{Arc, Mutex};

use anyhow::{self, Error};
use bincode::{self, deserialize_from, serialize_into};
use flate2::read::MultiGzDecoder;
use serde::{de::DeserializeOwned, Serialize};

use bio::io::{fasta, fastq};
use debruijn::dna_string::DnaString;
use log::info;

use crate::config::FastaFormat;

pub fn write_obj<T: Serialize, P: AsRef<Path> + Debug>(
    g: &T,
    filename: P,
) -> Result<(), bincode::Error> {
    let f = match File::create(&filename) {
        Err(err) => panic!("couldn't create file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut writer = BufWriter::new(f);
    serialize_into(&mut writer, &g)
}

pub fn read_obj<T: DeserializeOwned, P: AsRef<Path> + Debug>(
    filename: P,
) -> Result<T, bincode::Error> {
    let f = match File::open(&filename) {
        Err(err) => panic!("couldn't open file {:?}: {}", filename, err),
        Ok(f) => f,
    };
    let mut reader = BufReader::new(f);
    deserialize_from(&mut reader)
}

/// Open a (possibly gzipped) file into a BufReader.
fn _open_with_gz<P: AsRef<Path>>(p: P) -> Result<Box<dyn BufRead>, Error> {
    let r = File::open(p.as_ref())?;

    if p.as_ref().extension().unwrap() == "gz" {
        let gz = MultiGzDecoder::new(r);
        let buf_reader = BufReader::with_capacity(32 * 1024, gz);
        Ok(Box::new(buf_reader))
    } else {
        let buf_reader = BufReader::with_capacity(32 * 1024, r);
        Ok(Box::new(buf_reader))
    }
}

type ReadTranscriptResults = (Vec<DnaString>, Vec<String>, HashMap<String, String>);

pub fn read_transcripts(
    reader: fasta::Reader<BufReader<File>>,
) -> Result<ReadTranscriptResults, Error> {
    let mut seqs = Vec::new();
    let mut transcript_counter = 0;
    let mut tx_ids = Vec::new();
    let mut tx_to_gene_map = HashMap::new();
    let mut fasta_format = FastaFormat::Unknown;

    info!("Reading transcripts from Fasta file");
    for result in reader.records() {
        // obtain record or fail with error
        let record = result?;

        // Sequence
        let dna_string = DnaString::from_acgt_bytes_hashn(record.seq(), record.id().as_bytes());
        seqs.push(dna_string);

        if let FastaFormat::Unknown = fasta_format {
            fasta_format = detect_fasta_format(&record)?;
        }

        let (tx_id, gene_id) = extract_tx_gene_id(&record, &fasta_format);

        tx_ids.push(tx_id.to_string());
        tx_to_gene_map.insert(tx_id.to_string(), gene_id.to_string());

        transcript_counter += 1;
    }

    info!(
        "Done reading the Fasta file; Found {} sequences",
        transcript_counter
    );

    Ok((seqs, tx_ids, tx_to_gene_map))
}

fn detect_fasta_format(record: &fasta::Record) -> Result<FastaFormat, Error> {
    let id_tokens = record.id().split('|');
    if id_tokens.count() == 9 {
        return Ok(FastaFormat::Gencode);
    }

    let mut desc_tokens = record.desc().unwrap().split(' ');
    if let Some(desc_token) = desc_tokens.next() {
        let mut gene_tokens = desc_token.split('=');
        if let Some(gene_token) = gene_tokens.next() {
            if gene_token == "gene" && gene_tokens.count() == 1 {
                return Ok(FastaFormat::Gffread);
            }
        } else if desc_tokens.count() == 4 {
            return Ok(FastaFormat::Ensembl);
        }
    }
    anyhow::bail!("Failed to detect FASTA header format.")
}

fn extract_tx_gene_id<'a>(
    record: &'a fasta::Record,
    fasta_format: &FastaFormat,
) -> (&'a str, &'a str) {
    match *fasta_format {
        FastaFormat::Gencode => {
            let mut id_tokens = record.id().split('|');
            let tx_id = id_tokens.next().unwrap();
            let gene_id = id_tokens.next().unwrap();
            // (human readable name)
            // let gene_name = id_tokens[5].to_string();
            (tx_id, gene_id)
        }
        FastaFormat::Ensembl => {
            let tx_id = record.id();
            let mut desc_tokens = record.desc().unwrap().split(' ');
            let gene_id = desc_tokens.nth(2).unwrap().split(':').nth(1).unwrap();
            (tx_id, gene_id)
        }
        FastaFormat::Gffread => {
            let mut id_tokens = record.id().split(' ');
            let tx_id = id_tokens.next().unwrap();
            let mut desc_tokens = record.desc().unwrap().split(' ');
            let mut gene_tokens = desc_tokens.next().unwrap().split('=');
            let gene_id = gene_tokens.nth(1).unwrap();
            (tx_id, gene_id)
        }
        FastaFormat::Unknown => {
            panic!("fasta_format was uninitialized");
        }
    }
}

pub(crate) fn get_next_record<R: io::BufRead>(
    reader: &Arc<Mutex<fastq::Records<R>>>,
) -> Option<Result<fastq::Record, bio::io::fastq::Error>> {
    let mut lock = reader.lock().unwrap();
    lock.next()
}

pub(crate) fn open_file<P: AsRef<Path>>(filename: &str, outdir: P) -> Result<File, Error> {
    let out_fn = outdir.as_ref().join(filename);
    let outfile = File::create(&out_fn)?;
    Ok(outfile)
}
