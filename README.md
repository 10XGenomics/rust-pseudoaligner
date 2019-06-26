# Rust Pseudoaligner

A work-in-progress tool for pseudo-alignment of RNA-seq reads to transcriptome references. This project aims to create a very high performance pseudo-alignment tool suitable for single-cell RNA-seq data, and easily usable as a component of larger pipelines. We build on the crates includng [debruijn](https://github.com/10XGenomics/rust-debruijn) and [boomphf](https://github.com/10XGenomics/rust-boomphf), and [bio](https://github.com/rust-bio/rust-bio)

This tool implements existing algorithms from the literature inlcuding:

Bray, Nicolas L., et al. "Near-optimal probabilistic RNA-seq quantification." Nature biotechnology 34.5 (2016): 525.
https://arxiv.org/ftp/arxiv/papers/1505/1505.02710.pdf

Srivastava, Avi, et al. "RapMap: a rapid, sensitive and accurate tool for mapping RNA-seq reads to transcriptomes." Bioinformatics 32.12 (2016): i192-i200.
https://www.biorxiv.org/content/biorxiv/early/2015/10/22/029652.full.pdf
	
Ntranos, Vasilis, et al. "Fast and accurate single-cell RNA-seq analysis by clustering of transcript-compatibility counts." Genome biology 17.1 (2016): 112.
https://www.biorxiv.org/content/biorxiv/early/2016/03/04/036863.full.pdf
	
Srivastava, Avi, et al. "Alevin: An integrated method for dscRNA-seq quantification." bioRxiv (2018): 335000
https://www.biorxiv.org/content/biorxiv/early/2018/06/01/335000.full.pdf

Limasset, Antoine, et al. "Fast and scalable minimal perfect hashing for massive key sets." arXiv preprint arXiv:1702.03154 (2017).
https://arxiv.org/pdf/1702.03154.pdf

Orenstein, Yaron, et al. "Designing small universal k-mer hitting sets for improved analysis of high-throughput sequencing." PLoS computational biology 13.10 (2017): e1005777.
https://journals.plos.org/ploscompbiol/article/file?id=10.1371/journal.pcbi.1005777&type=printable

Li, Yang. "MSPKmerCounter: a fast and memory efficient approach for k-mer counting." arXiv preprint arXiv:1505.06550 (2015).
https://arxiv.org/pdf/1505.06550.pdf
