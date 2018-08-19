# Rust Pseudoaligner

A work-in-progress tool for pseudo-alignment of RNA-seq reads to transcriptome references. This project aims to create a very high performance pseudo-alignment tool suitable for single-cell RNA-seq data, and easily usable as a component of larger pipelines. We build on the crates includng[debruijn](https://github.com/10XGenomics/rust-debruijn) and [boomphf](https://github.com/10XGenomics/rust-boomphf), and [bio](https://github.com/rust-bio/rust-bio)

This tool implements existing algorithms from the literature inlcuding:

Bray, Nicolas L., et al. "Near-optimal probabilistic RNA-seq quantification." Nature biotechnology 34.5 (2016): 525.

Srivastava, Avi, et al. "RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes." Bioinformatics 32.12 (2016): i192-i200.
	
Ntranos, Vasilis, et al. "Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts." Genome biology 17.1 (2016): 112.
	
Srivastava, Avi, et al. "Alevin: An integrated method for dscRNA-seq quantification." bioRxiv (2018): 335000

Limasset, Antoine, et al. "Fast and scalable minimal perfect hashing for massive key sets." arXiv preprint arXiv:1702.03154 (2017).

Orenstein, Yaron, et al. "Designing small universal k-mer hitting sets for improved analysis of high-throughput sequencing." PLoS computational biology 13.10 (2017): e1005777.

Li, Yang. "MSPKmerCounter: a fast and memory efficient approach for k-mer counting." arXiv preprint arXiv:1505.06550 (2015).