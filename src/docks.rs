// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

//! Utilities for reading universal hitting set Kmer libraries. Hitting sets calculated from:
//! "Designing small universal k-mer hitting sets for improved analysis of high-throughput sequencing", Yaron Orenstein et. al

use debruijn::dna_string::DnaString;
use debruijn::msp::MspInterval;
use debruijn::{Kmer, Vmer};

use std::fs::File;
use std::io::{BufRead, BufReader};

use config::{DocksUhs, Minimizer, DOCKS_FILE, K, L};

// Docks Universal hitting sets
pub fn read_uhs() -> DocksUhs {
    info!("Starting reading Docks' Universal Hitting Set");
    let input = match File::open(DOCKS_FILE) {
        Err(err) => panic!("coulnd't open file {:?}: {}", DOCKS_FILE, err),
        Ok(handle) => handle,
    };

    let buffered = BufReader::new(input);
    let mut universal_hitting_set: DocksUhs = DocksUhs::new();

    let mut kmer_id: u16 = 0;
    for line in buffered.lines() {
        universal_hitting_set.insert(line.expect("couln't read line in uhs"), kmer_id);
        kmer_id += 1;
    }

    info!("Using K = {}, L = {}", K, L);
    info!(
        "Read total {:?} kmers as Docks-UHS from {}",
        universal_hitting_set.len(),
        DOCKS_FILE
    );
    universal_hitting_set
}

pub fn generate_msps(seq: &DnaString, uhs: &DocksUhs) -> Vec<MspInterval> {
    // Can't partition strings shorter than k
    assert!(seq.len() >= L);
    assert!(K <= 8);
    assert!(seq.len() < 1 << 32);

    // lambda to get an index of a kmer in uhs
    let uhs_idx = |i: usize| {
        let pi = seq.get_kmer::<Minimizer>(i);
        match uhs.get(&pi.to_string()) {
            Some(minimizer) => minimizer,
            None => panic!("Can't find kmer in uhs for {:?}", &pi.to_string()),
        }
    };

    // lamnda to to get index of the minimum uhs
    let min_uhs_pos = |i: usize, j: usize| if uhs_idx(i) <= uhs_idx(j) { i } else { j };

    let in_uhs = |i: usize| {
        let pi = seq.get_kmer::<Minimizer>(i);

        // no need to check for reverse complement since
        // kmers search in DOCKS is NOT canonicalised
        uhs.contains_key(&pi.to_string())
    };

    // lambda on range [start, stop]
    let find_min_pos_in_range = |start, stop| {
        let mut pos = start;
        let mut min_pos: Option<usize> = None;

        while pos < stop + 1 {
            if in_uhs(pos) {
                match min_pos {
                    Some(old_pos) => {
                        min_pos = Some(min_uhs_pos(old_pos, pos));
                    }
                    None => {
                        min_pos = Some(pos);
                    }
                };
            }
            pos += 1;
        }

        match min_pos {
            None => panic!(
                "Can't find a minimizer for {:?}",
                seq.slice(start, stop - 1)
            ),
            Some(pos) => pos,
        }
    };

    let seq_len = seq.len();
    let mut min_positions = Vec::with_capacity(16);
    let mut min_pos = find_min_pos_in_range(0, L - K);
    min_positions.push((0, min_pos));

    for i in 1..(seq_len - L + 1) {
        if i > min_pos {
            min_pos = find_min_pos_in_range(i, i + L - K);
            min_positions.push((i, min_pos));
        } else if in_uhs(i + L - K) {
            let test_min = min_uhs_pos(min_pos, i + L - K);

            if test_min != min_pos {
                min_pos = test_min;
                min_positions.push((i, min_pos));
            }
        }
    }

    let mut slices = Vec::with_capacity(min_positions.len());

    // Generate the slices of the final string
    for p in 0..min_positions.len() - 1 {
        let (start_pos, min_pos) = min_positions[p];
        let (next_pos, _) = min_positions[p + 1];

        let interval = MspInterval::new(
            *uhs_idx(min_pos),
            start_pos as u32,
            (next_pos + L - 1 - start_pos) as u16,
        );
        slices.push(interval);
    }

    let (last_pos, min_pos) = min_positions[min_positions.len() - 1];
    let last_interval = MspInterval::new(
        *uhs_idx(min_pos),
        last_pos as u32,
        (seq_len - last_pos) as u16,
    );
    slices.push(last_interval);

    slices
}
