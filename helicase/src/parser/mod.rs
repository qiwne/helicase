//! Parser types and traits.

mod fasta;
mod fastq;
mod fastx;
mod traits;

pub use fasta::*;
pub use fastq::*;
pub use fastx::*;
pub use traits::*;

pub enum Format {
    Fasta,
    Fastq,
}

pub enum Event {
    Record(usize),
    DnaChunk(usize),
}
