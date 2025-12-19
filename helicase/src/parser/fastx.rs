use super::*;
use crate::config::{advanced::*, *};
use crate::dna_format::*;
use crate::input::*;

/// A wrapper for [`FastaParser`] / [`FastqParser`] detecting the format at runtime.
pub struct FastxParser<'a, const CONFIG: Config>(Box<dyn ParserIter + 'a>);

impl<'a, const CONFIG: Config, I: InputData<'a> + 'a> FromInputData<'a, I>
    for FastxParser<'a, CONFIG>
{
    fn from_input(input: I) -> Self {
        match input.first_byte() {
            b'>' => Self(Box::new(FastaParser::<CONFIG, I>::from_input(input))),
            b'@' => Self(Box::new(FastqParser::<CONFIG, I>::from_input(input))),
            _ => panic!("Unknow input format"),
        }
    }
}

impl<'a, const CONFIG: Config> Parser for FastxParser<'a, CONFIG> {
    #[inline(always)]
    fn format(&self) -> Format {
        self.0.format()
    }

    #[inline(always)]
    fn clear_record(&mut self) {
        self.0.clear_record();
    }

    #[inline(always)]
    fn clear_chunk(&mut self) {
        self.0.clear_chunk();
    }

    #[inline(always)]
    fn get_header(&self) -> &[u8] {
        assert!(flag_is_set(CONFIG, COMPUTE_HEADER));
        self.0.get_header()
    }

    #[inline(always)]
    fn get_header_owned(&mut self) -> Vec<u8> {
        assert!(flag_is_set(CONFIG, COMPUTE_HEADER));
        self.0.get_header_owned()
    }

    #[inline(always)]
    fn get_quality(&self) -> Option<&[u8]> {
        assert!(flag_is_set(CONFIG, COMPUTE_QUALITY));
        self.0.get_quality()
    }

    #[inline(always)]
    fn get_quality_owned(&mut self) -> Option<Vec<u8>> {
        assert!(flag_is_set(CONFIG, COMPUTE_QUALITY));
        self.0.get_quality_owned()
    }

    #[inline(always)]
    fn get_dna_string(&self) -> &[u8] {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_STRING));
        self.0.get_dna_string()
    }

    #[inline(always)]
    fn get_dna_string_owned(&mut self) -> Vec<u8> {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_STRING));
        self.0.get_dna_string_owned()
    }

    #[inline(always)]
    fn get_dna_columnar(&self) -> &ColumnarDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR));
        self.0.get_dna_columnar()
    }

    #[inline(always)]
    fn get_dna_columnar_owned(&mut self) -> ColumnarDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_COLUMNAR));
        self.0.get_dna_columnar_owned()
    }

    #[inline(always)]
    fn get_dna_packed(&self) -> &PackedDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_PACKED));
        self.0.get_dna_packed()
    }

    #[inline(always)]
    fn get_dna_packed_owned(&mut self) -> PackedDNA {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_PACKED));
        self.0.get_dna_packed_owned()
    }

    #[inline(always)]
    fn get_dna_len(&self) -> usize {
        assert!(flag_is_set(CONFIG, COMPUTE_DNA_LEN));
        self.0.get_dna_len()
    }
}

impl<'a, const CONFIG: Config> Iterator for FastxParser<'a, CONFIG> {
    type Item = Event;

    #[inline(always)]
    fn next(&mut self) -> Option<Self::Item> {
        self.0.next()
    }
}
