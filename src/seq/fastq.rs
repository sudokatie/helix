use crate::seq::encode_sequence;
use flate2::read::GzDecoder;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Sequence {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

pub struct FastqReader<R: BufRead> {
    reader: R,
}

impl<R: BufRead> FastqReader<R> {
    pub fn new(reader: R) -> Self {
        Self { reader }
    }

    pub fn next_record(&mut self) -> std::io::Result<Option<Sequence>> {
        let mut header = String::new();
        if self.reader.read_line(&mut header)? == 0 {
            return Ok(None);
        }

        if !header.starts_with('@') {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Expected @ header",
            ));
        }

        let id = header[1..]
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();

        let mut seq_line = String::new();
        self.reader.read_line(&mut seq_line)?;
        let seq = encode_sequence(seq_line.trim().as_bytes());

        let mut plus = String::new();
        self.reader.read_line(&mut plus)?;

        let mut qual_line = String::new();
        self.reader.read_line(&mut qual_line)?;
        let qual = qual_line.trim().as_bytes().to_vec();

        Ok(Some(Sequence { id, seq, qual }))
    }
}

impl<R: BufRead> Iterator for FastqReader<R> {
    type Item = std::io::Result<Sequence>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_record() {
            Ok(Some(seq)) => Some(Ok(seq)),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

pub fn open_fastq<P: AsRef<Path>>(
    path: P,
) -> std::io::Result<FastqReader<BufReader<Box<dyn Read + Send>>>> {
    let file = File::open(path.as_ref())?;
    let reader: Box<dyn Read + Send> =
        if path.as_ref().extension().is_some_and(|e| e == "gz") {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };
    Ok(FastqReader::new(BufReader::new(reader)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_parse_single_record() {
        let data = "@read1\nACGT\n+\nIIII\n";
        let mut reader = FastqReader::new(BufReader::new(Cursor::new(data)));
        let rec = reader.next_record().unwrap().unwrap();
        assert_eq!(rec.id, "read1");
        assert_eq!(rec.seq.len(), 4);
        assert_eq!(rec.qual.len(), 4);
    }

    #[test]
    fn test_parse_multiple_records() {
        let data = "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let mut reader = FastqReader::new(BufReader::new(Cursor::new(data)));

        let rec1 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec1.id, "read1");

        let rec2 = reader.next_record().unwrap().unwrap();
        assert_eq!(rec2.id, "read2");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_empty_file() {
        let data = "";
        let mut reader = FastqReader::new(BufReader::new(Cursor::new(data)));
        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_iterator() {
        let data = "@read1\nACGT\n+\nIIII\n@read2\nTGCA\n+\nJJJJ\n";
        let reader = FastqReader::new(BufReader::new(Cursor::new(data)));
        let records: Vec<_> = reader.map(|r| r.unwrap()).collect();
        assert_eq!(records.len(), 2);
    }
}
