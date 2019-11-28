use std::sync::Mutex;

/// Ingest pairs of (index, value) from multiple threads and set
/// vec[index] = value efficiently.
pub struct ScatterToVec<'a, T> {
    slices: Vec<Mutex<&'a mut [T]>>,
    chunk_bit_size: usize,
    max_buf_size: usize,
}

const CHUNK_BITS: usize = 20;
const BUF_ELEMENTS: usize = 16;

impl<'a, T> ScatterToVec<'a, T> {
    /// Create a new scatterer that permits efficiently writing (index, value)
    /// tuples into `data` from multiple threads.
    pub fn new(data: &'a mut [T]) -> ScatterToVec<'a, T> {
        let mut slices = Vec::new();
        let sz = 1 << CHUNK_BITS;

        let mut rest = data;

        while rest.len() > sz {
            let (l, r) = rest.split_at_mut(sz);
            slices.push(Mutex::new(l));
            rest = r;
        }

        slices.push(Mutex::new(rest));

        ScatterToVec {
            slices,
            chunk_bit_size: CHUNK_BITS,
            max_buf_size: BUF_ELEMENTS,
        }
    }

    /// Create a writer handle. Each thread that produces values
    /// should give given it's own handle to write values with.
    pub fn handle(&'a self) -> ScatterHandle<'a, T> {
        let mut bufs = Vec::with_capacity(self.slices.len());
        for _ in 0..self.slices.len() {
            bufs.push(vec![]);
        }

        ScatterHandle {
            max_buf_size: self.max_buf_size,
            scatter: self,
            bufs,
        }
    }
}

/// A handle to write (index, value) pairs
/// into the target slice.
pub struct ScatterHandle<'a, T> {
    max_buf_size: usize,
    scatter: &'a ScatterToVec<'a, T>,
    bufs: Vec<Vec<(usize, T)>>,
}

impl<'a, T> ScatterHandle<'a, T> {
    /// Set data[index] = value in the data slice.
    pub fn write(&mut self, index: usize, value: T) {
        let chunk = index >> self.scatter.chunk_bit_size;
        let buf = &mut self.bufs[chunk];
        buf.push((index, value));

        // If we've filled this buffer, write out the values
        if buf.len() == self.max_buf_size {
            self.flush_chunk(chunk);
        }
    }

    fn flush_chunk(&mut self, chunk: usize) {
        let buf = &mut self.bufs[chunk];
        let mut slice_to_write = self.scatter.slices[chunk].lock().unwrap();

        for (index, value) in buf.drain(..) {
            let mask = (1 << self.scatter.chunk_bit_size) - 1;
            let slice_pos = mask & index;
            slice_to_write[slice_pos] = value;
        }
    }
}

impl<'a, T> Drop for ScatterHandle<'a, T> {
    fn drop(&mut self) {
        for i in 0..self.bufs.len() {
            self.flush_chunk(i);
        }
    }
}
