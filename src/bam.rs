
use rust_htslib::bam::Read;
use rust_htslib::bam::{Record, Reader};
use failure::Error;
use debruijn::dna_string::DnaString;
use smallvec::SmallVec;
pub struct BamSeqReader {
    reader: Reader,
    tmp_record: Record,
}

impl BamSeqReader {
    pub fn new(reader: Reader) -> BamSeqReader {
        BamSeqReader {
            reader,
            tmp_record: Record::new(),
        }
    }
}

type Barcode = SmallVec<[u8; 24]>;
type Umi = SmallVec<[u8; 16]>;

pub const PROC_BC_SEQ_TAG: &'static [u8]     = b"CB";
pub const PROC_UMI_SEQ_TAG: &'static [u8]    = b"UB";

struct BamRead {
    sequence: DnaString,
}


impl Iterator for BamSeqReader {
    type Item = Result<(DnaString, Barcode, Umi), Error>;

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

            // Get original read sequence from record.

            //let s = self.tmp_record.seq_len();
            //let st = DnaString::new(s);

            // FIXME - make this faster!!
            let mut dna = DnaString::from_acgt_bytes_hashn(&self.tmp_record.seq().as_bytes(), self.tmp_record.qname());

            if self.tmp_record.is_reverse() {
                dna = dna.reverse();
            }

            let bc = self.tmp_record.aux(PROC_BC_SEQ_TAG).map(|x| Barcode::from_slice(x.string()));
            let umi = self.tmp_record.aux(PROC_UMI_SEQ_TAG).map(|x| Umi::from_slice(x.string()));

            if bc.is_some() && umi.is_some() {
                return Some(Ok((dna, bc.unwrap(), umi.unwrap())));
            }
        }
    }
}