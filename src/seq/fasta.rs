use crate::seq::encode_sequence;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Reference {
    pub name: String,
    pub seq: Vec<u8>,
}

impl Reference {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

pub struct FastaReader<R: BufRead> {
    reader: R,
    current_header: Option<String>,
}

impl<R: BufRead> FastaReader<R> {
    pub fn new(mut reader: R) -> Self {
        let mut first_line = String::new();
        let _ = reader.read_line(&mut first_line);
        let header = if first_line.starts_with('>') {
            Some(
                first_line[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
            )
        } else {
            None
        };
        Self {
            reader,
            current_header: header,
        }
    }

    pub fn next_record(&mut self) -> std::io::Result<Option<Reference>> {
        let name = match self.current_header.take() {
            Some(h) => h,
            None => return Ok(None),
        };

        let mut seq = Vec::new();
        let mut line = String::new();

        loop {
            line.clear();
            if self.reader.read_line(&mut line)? == 0 {
                break;
            }

            if line.starts_with('>') {
                self.current_header = Some(
                    line[1..]
                        .split_whitespace()
                        .next()
                        .unwrap_or("")
                        .to_string(),
                );
                break;
            }

            seq.extend(encode_sequence(line.trim().as_bytes()));
        }

        Ok(Some(Reference { name, seq }))
    }
}

impl<R: BufRead> Iterator for FastaReader<R> {
    type Item = std::io::Result<Reference>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(r)) => Some(Ok(r)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub fn open_fasta<P: AsRef<Path>>(path: P) -> std::io::Result<FastaReader<BufReader<File>>> {
    let file = File::open(path)?;
    Ok(FastaReader::new(BufReader::new(file)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_single_sequence() {
        let data = ">chr1\nACGT\n";
        let mut reader = FastaReader::new(BufReader::new(Cursor::new(data)));
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.name, "chr1");
        assert_eq!(rec.seq.len(), 4);
    }

    #[test]
    fn test_parse_wrapped_sequence() {
        let data = ">chr1\nACGT\nTGCA\n";
        let mut reader = FastaReader::new(BufReader::new(Cursor::new(data)));
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.name, "chr1");
        assert_eq!(rec.seq.len(), 8);
    }

    #[test]
    fn test_parse_multiple_sequences() {
        let data = ">chr1\nACGT\n>chr2\nTGCA\n";
        let mut reader = FastaReader::new(BufReader::new(Cursor::new(data)));

        let r1 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.name, "chr1");

        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r2.name, "chr2");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_header_with_description() {
        let data = ">chr1 human chromosome 1\nACGT\n";
        let mut reader = FastaReader::new(BufReader::new(Cursor::new(data)));
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.name, "chr1"); // Only ID, not description
    }

    #[test]
    fn test_iterator() {
        let data = ">chr1\nACGT\n>chr2\nTGCA\n";
        let reader = FastaReader::new(BufReader::new(Cursor::new(data)));
        let refs: Vec<_> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(refs.len(), 2);
    }
}
