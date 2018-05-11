extern crate debruijn;

use std;
use std::marker::PhantomData;

use debruijn::{compression};

// Extending trait CompressionSpec for compression
pub struct ScmapCompress<D> {
    d: PhantomData<D>,
}

impl<D> ScmapCompress<D> {
    pub fn new() -> ScmapCompress<D> {
        ScmapCompress {
            d: PhantomData,
        }
    }
}

impl<D: PartialEq> compression::CompressionSpec<D> for ScmapCompress<D>
where D: std::fmt::Debug
{
    fn reduce(&self, d: D, other: &D) -> D {
        if d != *other {
            panic!("{:?} != {:?}, Should not happen", d, *other);
        }
        d
    }

    fn join_test(&self, d1: &D, d2: &D) -> bool {
        if d1 == d2 { true } else { false }
    }
}
