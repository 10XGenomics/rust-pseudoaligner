// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

use debruijn::Kmer;
use failure::Error;
use itertools::Itertools;

use config::MAPPABILITY_COUNTS_LEN;
use pseudoaligner::Pseudoaligner;

// 1. Given graph, build a data structure of transcripts
//    - tx: tx_name, gene_name, 
// 2. For each de Bruijn graph node
//    - count = number of kmers (L - K + 1)
//    - transcript multiplicity = # of colors (size of equiv class)
//    - gene multiplicity = # of distinct genes
//    - add count, transcript multiplicity to tx_mappability
//    - add count, gene multiplicity to gene_mappability
// 3. Output results to tx_mappability.tsv and gene_mappability.tsc
//    - tx_mappability:
//      tx_name gene_name length kmer_count fraction_unique_tx fraction_unique_gene
// MappabilityRecord: tx_name, gene_name, tx_multiplicity: [usize], gene_multiplicity: [usize]
//    
// fn update_counts(Vec<Record>, kmer_count, ids)
//    fn update_counts(self, kmer_count, ids), Option(Gene_tx_map))
//      - (if gene we'll need to make a gene vector instead of color)
//    fn fraction_unique(self) -> f64

#[derive(Debug)]
pub struct MappabilityRecord {
    pub tx_name: String,
    pub gene_name: String,
    tx_multiplicity: [usize; MAPPABILITY_COUNTS_LEN],
    gene_multiplicity: [usize; MAPPABILITY_COUNTS_LEN]
}

impl MappabilityRecord {
    pub fn new(tx_name: &String, gene_name: &String) -> MappabilityRecord {
        MappabilityRecord {
            tx_name: tx_name.clone(),
            gene_name: gene_name.clone(),
            // tx_multiplicity[j] = # of kmers in this tx shared by j other transcripts
            tx_multiplicity: [0; MAPPABILITY_COUNTS_LEN],
            // gene_multiplicity[j] = # of kmers in the tx shared by j other genes
            gene_multiplicity: [0; MAPPABILITY_COUNTS_LEN]
        }
    }

    pub fn total_kmer_count(&self) -> usize {
        self.tx_multiplicity
            .iter()
            .sum()
    }

    pub fn add_tx_count(&mut self, count: usize, multiplicity: usize) {
        if multiplicity > MAPPABILITY_COUNTS_LEN {
            self.tx_multiplicity[MAPPABILITY_COUNTS_LEN - 1] += count
        } else {
            self.tx_multiplicity[multiplicity - 1] += count
        }
    }

    pub fn add_gene_count(&mut self, count: usize, multiplicity: usize) {
        if multiplicity > MAPPABILITY_COUNTS_LEN {
            self.gene_multiplicity[MAPPABILITY_COUNTS_LEN - 1] += count
        } else {
            self.gene_multiplicity[multiplicity - 1] += count
        }
    }
    
    pub fn fraction_unique_tx(&self) -> f64 {
        self.tx_multiplicity[0] as f64 / self.total_kmer_count() as f64
    }

    pub fn fraction_unique_gene(&self) -> f64 {
        self.gene_multiplicity[0] as f64 / self.total_kmer_count() as f64
    }
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

pub fn analyze_graph<K: Kmer>(index: &Pseudoaligner<K>) -> Result<(Vec<MappabilityRecord>), Error> {
    let mut records = Vec::new();

    // Make records
    for tx_name in index.tx_names.iter() {
        let gene_name = index.tx_gene_mapping.get(tx_name).unwrap();
        records.push(MappabilityRecord::new(&tx_name, &gene_name));
    }

    // Iterate through graph
    for node in index.dbg.iter_nodes() {
        let num_kmer = node.len() - K::k() + 1;

        let eq_class_idx = *node.data() as usize;
        let eq_class = &index.eq_classes[eq_class_idx];

        let num_tx = eq_class.len();

        let mut gene_names = Vec::new();
        for &tx_id in eq_class {
            let tx_name = &index.tx_names[tx_id as usize];
            let gene_name = index.tx_gene_mapping.get(tx_name);
            gene_names.push(gene_name.clone())
        }
        let unique_genes: Vec<_> = gene_names.iter().unique().collect();
        let num_genes = unique_genes.len();

        for &tx_id in eq_class {
            records[tx_id as usize].add_tx_count(num_kmer, num_tx);
            records[tx_id as usize].add_gene_count(num_kmer, num_genes);
        }
    }

    Ok(records)
}
