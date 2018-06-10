// ! Reading in the universal hitting sets
// Designing small universal k-mer hitting sets for improved analysis of high-throughput sequencing Yaron Orenstein et. al

use debruijn::msp::MspInterval;
use debruijn::{Vmer, Kmer, kmer};
use debruijn::dna_string::DnaString;

use std::fs::File;
use std::collections::HashMap;
use std::io::{BufReader, BufRead};

// Docks Universal hitting sets
pub type DocksUhs = HashMap<String, u16>;
pub type Minimizer = kmer::Kmer8;

const K: usize = 8;
pub const L: usize = 80;
const DOCKS_FILE: &'static str = "res_8_80_4_0.txt";

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
        universal_hitting_set.insert(line.expect("couln't read line in uhs"),
                                     kmer_id.clone());
        kmer_id += 1;
    }

    info!("Read total {:?} kmers as Docks-UHS", universal_hitting_set.len());
    universal_hitting_set
}

pub fn generate_msps( seq: &DnaString, uhs: &DocksUhs)
                  -> Vec<MspInterval> {
    // Can't partition strings shorter than k
    assert!(seq.len() >= L);
    assert!(K <= 8);
    assert!(seq.len() < 1<<32);

    // lambda to get an index of a kmer in uhs
    let uhs_idx = | i: usize | {
        let pi = seq.get_kmer::<kmer::Kmer8>(i);
        uhs.get(&pi.to_string()).expect("Can't find kmer in uhs")
    };

    // lamnda to to get index of the minimum uhs
    let min_uhs_pos = |i: usize, j: usize| if uhs_idx(i) <= uhs_idx(j) { i } else { j };

    let in_uhs = |i: usize| {
        let pi = seq.get_kmer::<Minimizer>(i);

        // no need to check for reverse complement since
        // kmers search in DOCKS is NOT canonicalised
        uhs.contains_key(&pi.to_string())
    };

    // lambda on range [start, stop) [NOTE: stop is non-inclusive]
    let find_min_pos_in_range = |start, stop| {
        let mut pos = start;
        let mut min_pos: Option<usize> = None;

        while pos < stop {
            if in_uhs(pos) {
                match min_pos {
                    Some(old_pos) => { min_pos = Some(min_uhs_pos(old_pos, pos)); },
                    None => { min_pos = Some(pos); },
                };
            }
            pos = pos + 1;
        }

        min_pos.expect("Can't find a minimizer")
    };

    let seq_len = seq.len();
    let mut min_positions = Vec::with_capacity(16);
    let mut min_pos = find_min_pos_in_range(0, L - K);
    min_positions.push((0, min_pos));

    for i in 1..(seq_len - L + 1) {
        if i > min_pos {
            min_pos = find_min_pos_in_range(i, i + L - K);
            min_positions.push((i, min_pos));
        } else {
            if in_uhs(i + L - K) {
                let test_min = min_uhs_pos(min_pos, i + L - K);

                if test_min != min_pos {
                    min_pos = test_min;
                    min_positions.push((i, min_pos));
                }
            }
        }
    }

    let mut slices = Vec::with_capacity(min_positions.len());

    // Generate the slices of the final string
    for p in 0..min_positions.len() - 1 {
        let (start_pos, min_pos) = min_positions[p];
        let (next_pos, _) = min_positions[p + 1];

        let interval = MspInterval::new(
            uhs_idx(min_pos).clone() as u16,
            start_pos as u32,
            (next_pos + L - 1 - start_pos) as u16
        );
        slices.push(interval);
    }

    let (last_pos, min_pos) = min_positions[min_positions.len() - 1];
    let last_interval = MspInterval::new(
        uhs_idx(min_pos).clone() as u16,
        last_pos as u32,
        (seq_len - last_pos) as u16
    );
    slices.push(last_interval);

    slices
}


//    let mut msps: Vec<(u16, Exts, DnaString)> = Vec::new();
//        let seq_slice = &seq.to_bytes()[..];
//        for msp in msp_parts {
//            let v = DnaString::from_bytes(&seq_slice[(msp.start())..(msp.start() + msp.len())]);
//            let exts = Exts::from_slice_bounds(&seq_slice, msp.start(), msp.len());
//            msps.push((msp.bucket(), exts, v));
