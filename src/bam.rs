use debruijn::dna_string::DnaString;
use failure::Error;
use rust_htslib::bam::Read;
use rust_htslib::bam::{IndexedReader, Record};

use regex::Regex;
use shardio::{ShardWriter, SortKey};
use std::borrow::Cow;
use std::collections::HashMap;
use std::path::Path;
use std::str::FromStr;

use log::debug;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;

use crate::config::KmerType;
use crate::locus::Locus;
use crate::pseudoaligner::intersect;
use crate::pseudoaligner::Pseudoaligner;

pub struct BamSeqReader {
    reader: IndexedReader,
    tmp_record: Record,
    gene_regex: Regex,
}

pub const HLA_FILTER: &'static str = "^HLA-.*";

impl BamSeqReader {
    pub fn new(reader: IndexedReader) -> BamSeqReader {
        BamSeqReader {
            reader,
            tmp_record: Record::new(),
            gene_regex: Regex::new(HLA_FILTER).unwrap(),
        }
    }

    pub fn fetch(&mut self, locus: &Locus) {
        let tid = self.reader.header().tid(locus.chrom.as_bytes()).unwrap();
        self.reader.fetch(tid, locus.start, locus.end).unwrap();
    }
}

pub type Barcode = SmallVec<[u8; 24]>;
pub type Umi = SmallVec<[u8; 16]>;
pub type EqClass = SmallVec<[u32; 4]>;

pub const GENE_TAG: &'static [u8] = b"GN";
pub const PROC_BC_SEQ_TAG: &'static [u8] = b"CB";
pub const PROC_UMI_SEQ_TAG: &'static [u8] = b"UB";

#[derive(Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BamCrRead {
    key: BamCrReadKey,
    sequence: DnaString,
}

#[derive(Clone, Debug, Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BamCrReadKey {
    barcode: Barcode,
    umi: Umi,
}

struct Key;
impl SortKey<BamCrRead> for Key {
    type Key = BamCrReadKey;
    fn sort_key(item: &BamCrRead) -> Cow<BamCrReadKey> {
        Cow::Borrowed(&item.key)
    }
}

#[derive(Serialize, Deserialize, Ord, PartialOrd, Eq, PartialEq)]
pub struct BusCount {
    barcode_id: u32,
    umi_id: u32,
    eq_class_id: u32,
}

#[derive(Serialize, Deserialize)]
pub struct EqClassDb {
    barcodes: HashMap<Barcode, u32>,
    umis: HashMap<Umi, u32>,
    eq_classes: HashMap<EqClass, u32>,
    counts: Vec<BusCount>,
}

impl<'a> EqClassDb {
    pub fn new(_nitems: usize) -> EqClassDb {
        EqClassDb {
            barcodes: HashMap::new(),
            umis: HashMap::new(),
            eq_classes: HashMap::new(),
            counts: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.counts.len()
    }

    pub fn count(&mut self, barcode: &Barcode, umi: &Umi, eq_class: &EqClass) {
        let barcode_id = match self.barcodes.get(barcode).cloned() {
            Some(id) => id,
            None => {
                let id = self.barcodes.len() as u32;
                self.barcodes.insert(barcode.clone(), id);
                id
            }
        };

        let umi_id = match self.umis.get(umi).cloned() {
            Some(id) => id,
            None => {
                let id = self.umis.len() as u32;
                self.umis.insert(umi.clone(), id);
                id
            }
        };

        let eq_class_id = match self.eq_classes.get(eq_class).cloned() {
            Some(id) => id,
            None => {
                let id = self.eq_classes.len() as u32;
                self.eq_classes.insert(eq_class.clone(), id);
                id
            }
        };

        let count = BusCount {
            barcode_id,
            umi_id,
            eq_class_id,
        };

        self.counts.push(count);
    }

    /*
    pub fn get(&self, i: usize) -> (&Barcode, &Umi, &EqClass) {

        let ids = &self.counts[i];

        let bc = self.barcodes.get_by_right(&ids.barcode_id).unwrap();
        let umi = self.umis.get_by_right(&ids.umi_id).unwrap();
        let eq_class = self.eq_classes.get_by_right(&ids.eq_class_id).unwrap();

        (bc, umi, eq_class)
    }
    */

    pub fn sort(&mut self) {
        self.counts.sort();
    }

    pub fn eq_class_counts(&mut self, nitems: usize) -> crate::em::EqClassCounts {
        let mut rev_map = HashMap::<u32, &EqClass>::new();
        let mut counts: HashMap<EqClass, u32> = HashMap::new();

        self.sort();

        for (cls, id) in &self.eq_classes {
            rev_map.insert(*id, cls);
        }

        let mut uniq = 0;
        let mut total_reads = 0;
        let mut buf = Vec::new();

        use itertools::Itertools;
        for (_, mut hits) in &self.counts.iter().group_by(|c| (c.barcode_id, c.umi_id)) {
            let c = hits.next().unwrap();
            buf.clear();
            buf.extend(rev_map[&c.eq_class_id]);
            total_reads += 1;

            for c in hits {
                intersect(&mut buf, rev_map[&c.eq_class_id]);
                total_reads += 1;
            }

            let eqclass = EqClass::from_slice(&buf);
            let count = counts.entry(eqclass).or_default();
            *count += 1;

            uniq += 1;
        }

        println!("mean umis/read: {}", (uniq as f64) / (total_reads as f64));

        let empty = EqClass::new();
        println!(
            "empty eq_class counts: {:?} of {}",
            counts.get(&empty),
            self.counts.len()
        );

        crate::em::EqClassCounts {
            nitems: nitems,
            counts,
        }
    }
}

impl Iterator for BamSeqReader {
    type Item = Result<BamCrRead, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let r = self.reader.read(&mut self.tmp_record);

            match r {
                Ok(true) => (),
                Ok(false) => return None,
                Err(e) => return Some(Err(e.into())),
            }

            if self.tmp_record.is_secondary() || self.tmp_record.is_supplementary() {
                continue;
            }

            // Use reads that have no GN tag, or a GN tag that matches "HLA-*"
            let gene_filter = match self.tmp_record.aux(GENE_TAG) {
                Some(gn_aux) => {
                    let gn_bytes = gn_aux.string();

                    let gn_iter = gn_bytes.split(|x| *x == b';');
                    let mut result = false;
                    for gn in gn_iter {
                        let gn_str = std::str::from_utf8(gn).unwrap();

                        if self.gene_regex.is_match(gn_str) {
                            result = true;
                            break;
                        }
                    }
                    result
                }
                None => true,
            };

            if !gene_filter {
                continue;
            };

            // Get original read sequence from record.
            let mut sequence = DnaString::from_acgt_bytes_hashn(
                &self.tmp_record.seq().as_bytes(),
                self.tmp_record.qname(),
            );

            if !self.tmp_record.is_reverse() {
                sequence = sequence.reverse();
            }

            let barcode = match self
                .tmp_record
                .aux(PROC_BC_SEQ_TAG)
                .map(|x| Barcode::from_slice(x.string()))
            {
                Some(bc) => bc,
                None => continue,
            };

            let umi = match self
                .tmp_record
                .aux(PROC_UMI_SEQ_TAG)
                .map(|x| Umi::from_slice(x.string()))
            {
                Some(bc) => bc,
                None => continue,
            };

            return Some(Ok(BamCrRead {
                sequence,
                key: BamCrReadKey { barcode, umi },
            }));
        }
    }
}

pub fn map_bam(
    bam: impl AsRef<Path>,
    align: Pseudoaligner<KmerType>,
    locus_string: &Option<String>,
    outs: &Path,
) -> Result<(), Error> {
    let mut rdr = IndexedReader::from_path(bam)?;
    rdr.set_threads(4).unwrap();

    let mut itr = BamSeqReader::new(rdr);

    match locus_string {
        Some(l) => {
            let locus = Locus::from_str(l)?;
            println!("Locus: {:?}", locus);
            itr.fetch(&locus);
        }
        _ => (),
    }

    let mut hits_file = outs.to_path_buf();
    hits_file.set_extension("counts.bin");

    let mut reads_file = outs.to_path_buf();
    reads_file.set_extension("shards.bin");

    let _sw = ShardWriter::<BamCrRead, Key>::new(&reads_file, 16, 64, 1 << 20);

    let mut rid = 0;
    let mut some_aln = 0;
    let mut long_aln = 0;

    let mut eq_counts = EqClassDb::new(align.tx_names.len());
    let mut nodes = Vec::new();
    let mut eq_class = Vec::new();

    for _rec in itr {
        rid += 1;
        let rec = _rec?;

        let aln = align.map_read_to_nodes(&rec.sequence, &mut nodes);

        if let Some(cov) = aln {
            some_aln += 1;

            align.nodes_to_eq_class(&mut nodes, &mut eq_class);

            if cov > 50 {
                long_aln += 1;

                eq_counts.count(
                    &rec.key.barcode,
                    &rec.key.umi,
                    &EqClass::from_slice(&eq_class),
                );

                /*
                for n in &nodes {
                    let bc = std::str::from_utf8(&rec.barcode).unwrap();
                    let umi = std::str::from_utf8(&rec.umi).unwrap();
                    writeln!(hits, "{}\t{}\t{}\t{}\t{}", bc, umi, rid, n, cov)?;
                }
                */
            }
        }

        if rid % 100000 == 0 {
            println!(
                "analyzed {} reads. Mapped {}, long {}",
                rid, some_aln, long_aln
            );
        }
    }

    println!(
        "analyzed {} reads. Mapped {}, long {}",
        rid, some_aln, long_aln
    );

    /*
    for node in align.dbg.iter_nodes() {
        let eq_class_id = node.data();
        use debruijn::Mer;
        let seq_len = node.sequence().len();
        let eq_class = &align.eq_classes[*eq_class_id as usize];

        for tx_id in eq_class.iter() {
            writeln!(idx, "{}\t{}\t{}\t{}\t{}", node.node_id, eq_class.len(), seq_len, tx_id, align.tx_names[*tx_id as usize]);
        }
    }
    */

    crate::utils::write_obj(&eq_counts, &hits_file)?;
    Ok(())
}

fn filter_rec(rec: &Record, locus: &Locus) -> bool {
    if rec.mapq() < 30 {
        debug!(
            "skipping read {} due to low mapping quality",
            String::from_utf8(rec.qname().to_vec()).unwrap()
        );

        false
    } else if useful_alignment(locus, &rec).unwrap() == false {
        debug!(
            "skipping read {} due to not being useful",
            String::from_utf8(rec.qname().to_vec()).unwrap()
        );
        false
    } else {
        true
    }
}

pub fn useful_alignment(locus: &Locus, rec: &Record) -> Result<bool, Error> {
    // filter alignments to ensure that they truly overlap the region of interest
    // for now, overlap will be defined as having an aligned base anywhere in the locus
    let cigar = rec.cigar();
    for i in locus.start..=locus.end {
        // Don't include soft-clips but do include deletions
        let t = cigar.read_pos(i, false, true)?;
        if t.is_some() {
            return Ok(true);
        }
    }
    Ok(false)
}
