
use rust_htslib::bam::Read;
use rust_htslib::bam::{Record, IndexedReader};
use failure::Error;
use debruijn::dna_string::DnaString;
use pseudoaligner::Pseudoaligner;
use config::KmerType;
use std::path::Path;

use smallvec::SmallVec;

use crate::locus::Locus;

pub struct BamSeqReader {
    reader: IndexedReader,
    tmp_record: Record,

}


impl BamSeqReader  {
    pub fn new(reader: IndexedReader) -> BamSeqReader {
        BamSeqReader {
            reader,
            tmp_record: Record::new(),
        }
    }

    pub fn fetch(&mut self, locus: &Locus) {
        let tid = self.reader.header().tid(locus.chrom.as_bytes()).unwrap();
        self.reader.fetch(tid, locus.start, locus.end);
    }
}



type Barcode = SmallVec<[u8; 24]>;
type Umi = SmallVec<[u8; 16]>;

pub const PROC_BC_SEQ_TAG: &'static [u8]     = b"CB";
pub const PROC_UMI_SEQ_TAG: &'static [u8]    = b"UB";

#[derive(Debug)]
pub struct BamCrRead {
    sequence: DnaString,
    barcode: Barcode,
    umi: Umi,
}



impl Iterator for BamSeqReader {
    type Item = Result<BamCrRead, Error>;

    fn next(&mut self) -> Option<Self::Item> {

        loop {
            let r = self.reader.read(&mut self.tmp_record);

            match r {
                Err(e) => {
                    if e.is_eof() {
                        return None
                    } else {
                        return Some(Err(e.into()))
                    }
                },
                _ => (),
            };

            if self.tmp_record.is_secondary() || self.tmp_record.is_supplementary() {
                continue;
            }

            // Get original read sequence from record.

            //let s = self.tmp_record.seq_len();
            //let st = DnaString::new(s);

            // FIXME - make this faster!!
            let mut sequence = DnaString::from_acgt_bytes_hashn(&self.tmp_record.seq().as_bytes(), self.tmp_record.qname());

            if self.tmp_record.is_reverse() {
                sequence = sequence.reverse();
            }

            let barcode = match self.tmp_record.aux(PROC_BC_SEQ_TAG).map(|x| Barcode::from_slice(x.string())) {
                Some(bc) => bc,
                None => continue,
            };

            let umi = match self.tmp_record.aux(PROC_UMI_SEQ_TAG).map(|x| Umi::from_slice(x.string())) {
                Some(bc) => bc,
                None => continue,
            };


            return Some(Ok(BamCrRead { sequence, barcode, umi }));
        }
    }
}

use std::io::BufWriter;
use std::fs::File;
use std::io::Write;
use std::str::FromStr;

pub fn map_bam(bam: impl AsRef<Path>, align: Pseudoaligner<KmerType>, locus_string: &str, outs: &Path) -> Result<(), Error> {

    let locus = Locus::from_str(locus_string)?;
    let rdr = IndexedReader::from_path(bam)?;
    println!("Locus: {:?}", locus);

    let mut itr = BamSeqReader::new(rdr);
    itr.fetch(&locus);

    let mut hits_tsv = outs.to_path_buf();
    hits_tsv.set_extension("blah1.tsv");

    let mut idx_tsv = outs.to_path_buf();
    idx_tsv.set_extension("blah2.tsv");

    let mut hits = BufWriter::new(File::create(hits_tsv)?);
    let mut idx = BufWriter::new(File::create(idx_tsv)?);

    let mut rid = 0;

    for _rec in itr {
        rid += 1;
        let rec = _rec?;

        let aln = align.map_read_to_nodes(&rec.sequence);

        println!("rec: {:?}", rec);
        println!("res: {:?}", aln);

        if let Some((nodes, cov)) = aln {
            for n in nodes {
                let bc = std::str::from_utf8(&rec.barcode).unwrap();
                let umi = std::str::from_utf8(&rec.umi).unwrap();
                writeln!(hits, "{}\t{}\t{}\t{}\t{}", bc, umi, rid, n, cov)?;
            }
        }
    }

    for node in align.dbg.iter_nodes() {
        let eq_class_id = node.data();
        use debruijn::Mer;
        let seq_len = node.sequence().len();
        let eq_class = &align.eq_classes[*eq_class_id as usize];

        for tx_id in eq_class.iter() {
            writeln!(idx, "{}\t{}\t{}\t{}\t{}", node.node_id, eq_class.len(), seq_len, tx_id, align.tx_names[*tx_id as usize]);
        }
    }

    Ok(())
}


fn filter_rec(rec: &Record, locus: &Locus) -> bool {

    if rec.mapq() < 30 {
        debug!("skipping read {} due to low mapping quality", 
                String::from_utf8(rec.qname().to_vec()).unwrap());

        false
    }
    else if useful_alignment(locus, &rec).unwrap() == false {
        debug!("skipping read {} due to not being useful", 
                String::from_utf8(rec.qname().to_vec()).unwrap());
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
                return Ok(true)
            }
        }
        Ok(false)
}
