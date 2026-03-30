use std::io::{Result, Write};

// SAM flags
pub const FLAG_PAIRED: u16 = 0x1;
pub const FLAG_PROPER_PAIR: u16 = 0x2;
pub const FLAG_UNMAPPED: u16 = 0x4;
pub const FLAG_MATE_UNMAPPED: u16 = 0x8;
pub const FLAG_REVERSE: u16 = 0x10;
pub const FLAG_MATE_REVERSE: u16 = 0x20;
pub const FLAG_READ1: u16 = 0x40;
pub const FLAG_READ2: u16 = 0x80;
pub const FLAG_SECONDARY: u16 = 0x100;
pub const FLAG_SUPPLEMENTARY: u16 = 0x800;

/// SAM file header
#[derive(Debug, Clone)]
pub struct SamHeader {
    /// SAM format version
    pub version: String,
    /// Reference sequences: (name, length)
    pub references: Vec<(String, usize)>,
    /// Optional program info
    pub program: Option<String>,
}

impl Default for SamHeader {
    fn default() -> Self {
        Self {
            version: "1.6".to_string(),
            references: Vec::new(),
            program: None,
        }
    }
}

impl SamHeader {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_reference(&mut self, name: String, length: usize) {
        self.references.push((name, length));
    }
}

/// A single SAM alignment record
#[derive(Debug, Clone)]
pub struct SamRecord {
    /// Query name
    pub qname: String,
    /// Bitwise flags
    pub flag: u16,
    /// Reference name
    pub rname: String,
    /// 1-based position
    pub pos: u32,
    /// Mapping quality
    pub mapq: u8,
    /// CIGAR string
    pub cigar: String,
    /// Mate reference name
    pub rnext: String,
    /// Mate position
    pub pnext: u32,
    /// Template length
    pub tlen: i32,
    /// Sequence
    pub seq: String,
    /// Quality string
    pub qual: String,
    /// Optional tags
    pub tags: Vec<(String, String)>,
}

impl SamRecord {
    pub fn new_unmapped(qname: String, seq: String, qual: String) -> Self {
        Self {
            qname,
            flag: FLAG_UNMAPPED,
            rname: "*".to_string(),
            pos: 0,
            mapq: 0,
            cigar: "*".to_string(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq,
            qual,
            tags: Vec::new(),
        }
    }

    pub fn add_tag(&mut self, tag: &str, value: &str) {
        self.tags.push((tag.to_string(), value.to_string()));
    }
}

/// Writer for SAM format
pub struct SamWriter<W: Write> {
    writer: W,
}

impl<W: Write> SamWriter<W> {
    pub fn new(writer: W) -> Self {
        Self { writer }
    }

    /// Write the SAM header
    pub fn write_header(&mut self, header: &SamHeader) -> Result<()> {
        // HD line
        writeln!(self.writer, "@HD\tVN:{}\tSO:unsorted", header.version)?;

        // SQ lines
        for (name, length) in &header.references {
            writeln!(self.writer, "@SQ\tSN:{}\tLN:{}", name, length)?;
        }

        // PG line
        if let Some(ref program) = header.program {
            writeln!(self.writer, "@PG\tID:helix\tPN:{}", program)?;
        }

        Ok(())
    }

    /// Write a single alignment record
    pub fn write_record(&mut self, record: &SamRecord) -> Result<()> {
        write!(
            self.writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            record.qname,
            record.flag,
            record.rname,
            record.pos,
            record.mapq,
            record.cigar,
            record.rnext,
            record.pnext,
            record.tlen,
            record.seq,
            record.qual
        )?;

        // Write optional tags
        for (tag, value) in &record.tags {
            write!(self.writer, "\t{}", tag)?;
            write!(self.writer, ":{}", value)?;
        }

        writeln!(self.writer)?;
        Ok(())
    }

    /// Flush the writer
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn test_write_header() {
        let mut header = SamHeader::new();
        header.add_reference("chr1".to_string(), 1000);
        header.add_reference("chr2".to_string(), 2000);

        let mut buf = Vec::new();
        let mut writer = SamWriter::new(Cursor::new(&mut buf));
        writer.write_header(&header).unwrap();

        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("@HD\tVN:1.6"));
        assert!(output.contains("@SQ\tSN:chr1\tLN:1000"));
        assert!(output.contains("@SQ\tSN:chr2\tLN:2000"));
    }

    #[test]
    fn test_write_record() {
        let record = SamRecord {
            qname: "read1".to_string(),
            flag: 0,
            rname: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "50M".to_string(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
            tags: vec![],
        };

        let mut buf = Vec::new();
        let mut writer = SamWriter::new(Cursor::new(&mut buf));
        writer.write_record(&record).unwrap();

        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("read1\t0\tchr1\t100\t60\t50M\t*\t0\t0\tACGT\tIIII"));
    }

    #[test]
    fn test_unmapped_record() {
        let record = SamRecord::new_unmapped(
            "read1".to_string(),
            "ACGT".to_string(),
            "IIII".to_string(),
        );

        assert_eq!(record.flag, FLAG_UNMAPPED);
        assert_eq!(record.rname, "*");
        assert_eq!(record.cigar, "*");
    }

    #[test]
    fn test_record_with_tags() {
        let mut record = SamRecord {
            qname: "read1".to_string(),
            flag: 0,
            rname: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "50M".to_string(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq: "ACGT".to_string(),
            qual: "IIII".to_string(),
            tags: vec![],
        };
        record.add_tag("NM", "i:2");
        record.add_tag("AS", "i:95");

        let mut buf = Vec::new();
        let mut writer = SamWriter::new(Cursor::new(&mut buf));
        writer.write_record(&record).unwrap();

        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("NM:i:2"));
        assert!(output.contains("AS:i:95"));
    }

    #[test]
    fn test_flags() {
        assert_eq!(FLAG_PAIRED, 0x1);
        assert_eq!(FLAG_UNMAPPED, 0x4);
        assert_eq!(FLAG_REVERSE, 0x10);
        assert_eq!(FLAG_SECONDARY, 0x100);
    }

    #[test]
    fn test_complete_sam_file() {
        let mut header = SamHeader::new();
        header.add_reference("chr1".to_string(), 1000);

        let record = SamRecord {
            qname: "read1".to_string(),
            flag: 0,
            rname: "chr1".to_string(),
            pos: 100,
            mapq: 60,
            cigar: "50M".to_string(),
            rnext: "*".to_string(),
            pnext: 0,
            tlen: 0,
            seq: "ACGTACGTACGT".to_string(),
            qual: "IIIIIIIIIIII".to_string(),
            tags: vec![],
        };

        let mut buf = Vec::new();
        {
            let mut writer = SamWriter::new(Cursor::new(&mut buf));
            writer.write_header(&header).unwrap();
            writer.write_record(&record).unwrap();
        }

        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<_> = output.lines().collect();
        
        // Should have header and record
        assert!(lines.len() >= 2);
        assert!(lines[0].starts_with("@HD"));
        assert!(lines[1].starts_with("@SQ"));
        assert!(lines[2].starts_with("read1"));
    }
}
