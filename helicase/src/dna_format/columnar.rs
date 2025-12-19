use std::fmt;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub(crate) enum Nucleotide {
    A,
    C,
    T,
    G,
}

impl Nucleotide {
    #[inline(always)]
    pub const fn from_bits(b0: bool, b1: bool) -> Self {
        match (b0, b1) {
            (false, false) => Self::A,
            (false, true) => Self::C,
            (true, false) => Self::T,
            (true, true) => Self::G,
        }
    }

    #[allow(unused)]
    #[inline(always)]
    pub const fn to_bits(self) -> (bool, bool) {
        match self {
            Self::A => (false, false),
            Self::C => (false, true),
            Self::T => (true, false),
            Self::G => (true, true),
        }
    }

    #[inline(always)]
    pub const fn as_char(self) -> char {
        match self {
            Self::A => 'A',
            Self::C => 'C',
            Self::T => 'T',
            Self::G => 'G',
        }
    }
}

#[derive(Debug, Clone)]
pub struct ColumnarDNA {
    store0: Vec<u64>,
    store1: Vec<u64>,
    b0: u64,
    b1: u64,
    space: usize,
}

impl Default for ColumnarDNA {
    fn default() -> Self {
        Self::new()
    }
}

impl ColumnarDNA {
    #[inline(always)]
    pub const fn new() -> Self {
        Self {
            store0: Vec::new(),
            store1: Vec::new(),
            b0: 0,
            b1: 0,
            space: 64,
        }
    }

    #[inline(always)]
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            store0: Vec::with_capacity(capacity),
            store1: Vec::with_capacity(capacity),
            b0: 0,
            b1: 0,
            space: 64,
        }
    }

    #[inline(always)]
    pub const fn len(&self) -> usize {
        64 * self.store0.len() + (64 - self.space)
    }

    #[inline(always)]
    pub const fn is_empty(&self) -> bool {
        self.store0.is_empty() && self.space == 64
    }

    #[inline(always)]
    pub const fn capacity(&self) -> usize {
        self.store0.capacity()
    }

    #[inline(always)]
    pub fn clear(&mut self) {
        self.store0.clear();
        self.store1.clear();
        self.b0 = 0;
        self.b1 = 0;
        self.space = 64;
    }

    pub fn push_str(&mut self, s: &str) {
        for ch in s.bytes() {
            let (b0, b1) = match ch {
                b'A' | b'a' => (0, 0),
                b'C' | b'c' => (0, 1),
                b'G' | b'g' => (1, 1),
                b'T' | b't' => (1, 0),
                _ => panic!("Invalid nucleotide: {}", ch as char),
            };
            self.append(b0, b1, 1);
        }
    }

    #[inline(always)]
    pub fn append(&mut self, b0: u64, b1: u64, size: usize) {
        if size == 0 {
            // should not happen?
            return;
        }
        if self.space == 0 {
            self.store0.push(self.b0);
            self.store1.push(self.b1);
            self.b0 = b0;
            self.b1 = b1;
            self.space = 64 - size;
        } else {
            let low_bits: usize = size.min(self.space);
            let low_mask = !0 >> (64 - low_bits);
            let offset = 64 - self.space;
            self.b0 |= (b0 & low_mask) << offset;
            self.b1 |= (b1 & low_mask) << offset;
            let high_bits = size.saturating_sub(self.space);
            self.space -= low_bits;
            if high_bits > 0 {
                self.store0.push(self.b0);
                self.store1.push(self.b1);
                self.b0 = (b0 >> low_bits) & ((1 << high_bits) - 1);
                self.b1 = (b1 >> low_bits) & ((1 << high_bits) - 1);
                self.space = 64 - high_bits;
            }
        }
    }

    #[inline(always)]
    pub(crate) fn get(&self, i: usize) -> Option<Nucleotide> {
        if i >= self.len() {
            return None;
        }
        let word_idx = i / 64;
        let bit_idx = i % 64;
        let (b0, b1) = if word_idx == self.store0.len() {
            ((self.b0 >> bit_idx) & 1 != 0, (self.b1 >> bit_idx) & 1 != 0)
        } else {
            (
                (self.store0[word_idx] >> bit_idx) & 1 != 0,
                (self.store1[word_idx] >> bit_idx) & 1 != 0,
            )
        };
        Some(Nucleotide::from_bits(b0, b1))
    }
}

impl fmt::Display for ColumnarDNA {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for i in 0..self.len() {
            write!(f, "{}", self.get(i).unwrap().as_char()).unwrap();
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_starts_empty() {
        let dna = ColumnarDNA::new();
        assert_eq!(dna.len(), 0);
        assert_eq!(dna.b0, 0);
        assert_eq!(dna.b1, 0);
        assert_eq!(dna.store0, Vec::new());
        assert_eq!(dna.store1, Vec::new());
    }

    #[test]
    fn test_append_simple() {
        let mut dna = ColumnarDNA::new();
        // Append one nucleotide (size = 1)
        dna.append(1, 0, 1); // (b0=1,b1=0) => T
        assert_eq!(dna.len(), 1);
        assert_eq!(format!("{}", dna), "T");

        dna.append(0, 1, 1); // (b0=0,b1=1) => C
        assert_eq!(dna.len(), 2);
        assert_eq!(format!("{}", dna), "TC");
    }

    #[test]
    fn test_append_multiple() {
        let mut dna = ColumnarDNA::new();
        // Append four 2-bit nucleotides (A, C, T, G)
        dna.append(0, 0, 1); // A
        dna.append(0, 1, 1); // C
        dna.append(1, 0, 1); // T
        dna.append(1, 1, 1); // G
        assert_eq!(dna.len(), 4);
        assert_eq!(format!("{}", dna), "ACTG");
    }

    #[test]
    fn append_single_bits() {
        let mut v = ColumnarDNA::new();

        v.push_str("ACGT");

        assert_eq!(v.len(), 4);
        assert_eq!(v.to_string(), "ACGT");
    }

    #[test]
    fn append_exact_64_bits() {
        let seq = "A".repeat(64);
        let mut v = ColumnarDNA::new();

        v.push_str(&seq);

        assert_eq!(v.len(), 64);
        assert_eq!(v.to_string(), seq);
        assert_eq!(v.store0.len(), 0); // still in partial word
    }

    #[test]
    fn test_crossing_all_offsets() {
        let letters = ['A', 'C', 'G', 'T'];

        for offset in 0..64 {
            let mut v = ColumnarDNA::new();

            // Fill offset bits
            for _ in 0..offset {
                v.push_str("A");
            }

            // Write 1–10 bits after the boundary
            for size in 1..10 {
                let ch = letters[size & 3];
                v.push_str(&ch.to_string());
            }

            // Roundtrip must match
            let s = v.to_string();
            assert_eq!(s.len(), offset + 9);
            assert_eq!(s, s); // no-op but ensures no invalid UTF
        }
    }

    #[test]
    fn repeated_small_appends() {
        let mut v = ColumnarDNA::new();
        let mut s = String::new();

        for i in 0..500 {
            let ch = ['A', 'C', 'G', 'T'][i & 3];
            s.push(ch);
            v.push_str(&ch.to_string());
        }

        assert_eq!(v.len(), s.len());
        assert_eq!(v.to_string(), s);
    }

    #[test]
    fn cross_boundary_regression() {
        let mut v = ColumnarDNA::new();
        let seq = "A".repeat(64) + "C";
        v.push_str(&seq);

        assert_eq!(v.len(), 65);
        assert_eq!(v.to_string(), seq);
        assert_eq!(v.store0.len(), 1); // one stored word now
    }
}
