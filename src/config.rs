// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.

// transcriptome fasta header formats
pub enum FastaFormat {
    Unknown,
    Gencode,
    Ensembl,
    Gffread,
}

// main configs
pub const MEM_SIZE: usize = 1;
pub const MIN_KMERS: usize = 1;
pub const STRANDED: bool = true;
pub const REPORT_ALL_KMER: bool = false;
pub const READ_COVERAGE_THRESHOLD: usize = 32;
pub const LEFT_EXTEND_FRACTION: f64 = 0.2;
pub const DEFAULT_ALLOWED_MISMATCHES: usize = 2;

pub const U32_MAX: usize = u32::max_value() as usize;

// Transcriptome mappability
pub const MAPPABILITY_COUNTS_LEN: usize = 11;
