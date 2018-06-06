// ! Reading in the universal hitting sets
// Designing small universal k-mer hitting sets for improved analysis of high-throughput sequencing Yaron Orenstein et. al

use debruijn::{Exts, Vmer, Kmer, kmer, DnaSlice};
use debruijn::msp::MspInterval;
use debruijn::dna_string::DnaString;

use std::fs::File;
use std::cmp::min;
use std::io::{BufReader, BufRead};
use std::collections::HashMap;

// Docks Universal hitting sets
pub type DocksUhs = HashMap<String, u16>;
pub type Minimizer = kmer::Kmer8;

pub fn read_uhs(path: String) -> DocksUhs {
    info!("Starting reading Docks' Universal Hitting Set");
    let input = match File::open(&path) {
        Err(err) => panic!("coulnd't open file {:?}: {}", path, err),
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

fn generate_msps( l: usize, k: usize,
                  seq: &DnaString, uhs: &DocksUhs,
                  rc: bool)
                  -> Vec<MspInterval> {
    // Can't partition strings shorter than k
    assert!(seq.len() >= l);
    assert!(k <= 8);
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
    let mut min_pos = find_min_pos_in_range(0, l - k);
    min_positions.push((0, min_pos));

    for i in 1..(seq_len - l + 1) {
        if i > min_pos {
            min_pos = find_min_pos_in_range(i, i + l - k);
            min_positions.push((i, min_pos));
        } else {
            if in_uhs(i + l - k) {
                let test_min = min_uhs_pos(min_pos, i + l - k);

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
            (next_pos + l - 1 - start_pos) as u16
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


pub fn msp_sequence( seq: DnaString,
                     uhs: &DocksUhs)
                     -> Vec<(u16, Exts, DnaString)> {
    let l: usize = 40;
    let k: usize = 8;

    let mut msps: Vec<(u16, Exts, DnaString)> = Vec::new();
    if seq.len() >= l {
        // Check for the length of contigs broken between two `N`.

        let msp_parts = generate_msps( l, k , &seq,
                                       &uhs, true );

        let seq_slice = &seq.to_bytes()[..];
        for msp in msp_parts {
            let v = DnaString::from_bytes(&seq_slice[(msp.start())..(msp.start() + msp.len())]);
            let exts = Exts::from_slice_bounds(&seq_slice, msp.start(), msp.len());
            msps.push((msp.bucket(), exts, v));
        }
    }

    msps
}
