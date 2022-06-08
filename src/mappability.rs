// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use anyhow::Error;
use debruijn::Kmer;
use itertools::Itertools;
use std::io::Write;
use std::path::Path;

use crate::config::MAPPABILITY_COUNTS_LEN;
use crate::pseudoaligner::Pseudoaligner;
use crate::utils::open_file;

// 1. Given graph, build a data structure of transcripts
//    - tx: tx_name, gene_name,
// 2. For each de Bruijn graph node
//    - count = number of kmers (L - K + 1)
//    - transcript multiplicity = # of colors (size of equiv class)
//    - gene multiplicity = # of distinct genes
//    - add count, transcript multiplicity to tx_mappability
//    - add count, gene multiplicity to gene_mappability
// 3. Output results to tx_mappability.tsv and gene_mappability.tsv
//    - tx_mappability:
//      tx_name gene_name length kmer_count fraction_unique_tx fraction_unique_gene
// MappabilityRecord: tx_name, gene_name, tx_multiplicity: [usize], gene_multiplicity: [usize]
//
// fn update_counts(Vec<Record>, kmer_count, ids)
//    fn update_counts(self, kmer_count, ids), Option(Gene_tx_map))
//      - (if gene we'll need to make a gene vector instead of color)
//    fn fraction_unique(self) -> f64
const MAPPABILITY_HEADER_STRING: &str =
    "tx_name\tgene_name\ttx_kmer_count\tfrac_kmer_unique_tx\tfrac_kmer_unique_gene\n";

#[derive(Debug)]
pub struct MappabilityRecord<'a> {
    pub tx_name: &'a str,
    pub gene_name: &'a str,
    tx_multiplicity: [usize; MAPPABILITY_COUNTS_LEN],
    gene_multiplicity: [usize; MAPPABILITY_COUNTS_LEN],
}

impl MappabilityRecord<'_> {
    fn new<'a>(tx_name: &'a str, gene_name: &'a str) -> MappabilityRecord<'a> {
        MappabilityRecord {
            tx_name,
            gene_name,
            // tx_multiplicity[j] = # of kmers in this tx shared by j other transcripts
            tx_multiplicity: [0; MAPPABILITY_COUNTS_LEN],
            // gene_multiplicity[j] = # of kmers in the tx shared by j other genes
            gene_multiplicity: [0; MAPPABILITY_COUNTS_LEN],
        }
    }

    fn total_kmer_count(&self) -> usize {
        self.tx_multiplicity.iter().sum()
    }

    fn add_tx_count(&mut self, count: usize, multiplicity: usize) {
        if multiplicity > MAPPABILITY_COUNTS_LEN {
            self.tx_multiplicity[MAPPABILITY_COUNTS_LEN - 1] += count
        } else {
            self.tx_multiplicity[multiplicity - 1] += count
        }
    }

    fn add_gene_count(&mut self, count: usize, multiplicity: usize) {
        if multiplicity > MAPPABILITY_COUNTS_LEN {
            self.gene_multiplicity[MAPPABILITY_COUNTS_LEN - 1] += count
        } else {
            self.gene_multiplicity[multiplicity - 1] += count
        }
    }

    fn fraction_unique_tx(&self) -> f64 {
        self.tx_multiplicity[0] as f64 / self.total_kmer_count() as f64
    }

    fn fraction_unique_gene(&self) -> f64 {
        self.gene_multiplicity[0] as f64 / self.total_kmer_count() as f64
    }

    fn to_tsv(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}",
            self.tx_name,
            self.gene_name,
            self.total_kmer_count(),
            self.fraction_unique_tx(),
            self.fraction_unique_gene()
        )
    }
}

pub fn write_mappability_tsv<P: AsRef<Path>>(
    records: Vec<MappabilityRecord>,
    outdir: P,
) -> Result<(), Error> {
    let mut outfile = open_file("tx_mappability.tsv", outdir)?;

    outfile.write_all(MAPPABILITY_HEADER_STRING.as_bytes())?;

    for record in records {
        writeln!(outfile, "{}", record.to_tsv())?;
    }

    Ok(())
}

// pub fn update_counts(records: &mut Vec<MappabilityRecord>,
//                      kmer_count: usize,
//                      ids: Vec<usize>) {
//     let num_ids = ids.len();
//     for id in ids {
//         // add to total # kmers
//         records[id].add_count(kmer_count, 0);
//         // add to counts according to the size of this equiv class
//         records[id].add_count(kmer_count, num_ids)
//     }
// }

pub fn analyze_graph<K: Kmer>(
    index: &Pseudoaligner<K>,
) -> Result<Vec<MappabilityRecord<'_>>, Error> {
    // Make records
    let mut records = index
        .tx_names
        .iter()
        .map(|tx_name| MappabilityRecord::new(tx_name, index.tx_gene_mapping.get(tx_name).unwrap()))
        .collect::<Vec<_>>();

    // Iterate through graph
    for node in index.dbg.iter_nodes() {
        let num_kmer = node.len() - K::k() + 1;

        let eq_class_idx = *node.data() as usize;
        let eq_class = &index.eq_classes[eq_class_idx];

        let num_tx = eq_class.len();

        let num_genes = eq_class
            .iter()
            .map(|&tx_id| {
                let tx_name = &index.tx_names[tx_id as usize];
                index.tx_gene_mapping.get(tx_name)
            })
            .unique()
            .count();

        for &tx_id in eq_class {
            let record = &mut records[tx_id as usize];
            record.add_tx_count(num_kmer, num_tx);
            record.add_gene_count(num_kmer, num_genes);
        }
    }

    Ok(records)
}
