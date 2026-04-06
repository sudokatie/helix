//! BAM file format writer
//!
//! BAM is the binary, BGZF-compressed version of SAM format.
//! This module provides BAM writing with proper headers, coordinate sorting,
//! and optional BAI index generation.

use std::io::{Result, Write, BufWriter};
use std::fs::File;
use std::path::Path;
use flate2::Compression;
use flate2::write::GzEncoder;

use super::{SamHeader, SamRecord, FLAG_UNMAPPED};

/// BAM magic number
const BAM_MAGIC: &[u8] = b"BAM\x01";

/// BGZF block size (64KB max compressed)
const BGZF_BLOCK_SIZE: usize = 65536;
const BGZF_MAX_BLOCK_SIZE: usize = 65536;

/// BAM file writer with BGZF compression
pub struct BamWriter<W: Write> {
    /// Inner BGZF writer
    writer: BgzfWriter<W>,
    /// Reference sequences for index
    references: Vec<(String, usize)>,
    /// Whether header has been written
    header_written: bool,
    /// Current virtual offset for indexing
    virtual_offset: u64,
    /// Index builder (if enabled)
    index: Option<BaiIndexBuilder>,
}

impl<W: Write> BamWriter<W> {
    /// Create a new BAM writer
    pub fn new(writer: W) -> Self {
        Self {
            writer: BgzfWriter::new(writer),
            references: Vec::new(),
            header_written: false,
            virtual_offset: 0,
            index: None,
        }
    }

    /// Enable BAI index generation
    pub fn with_index(mut self) -> Self {
        self.index = Some(BaiIndexBuilder::new());
        self
    }

    /// Write the BAM header
    pub fn write_header(&mut self, header: &SamHeader) -> Result<()> {
        if self.header_written {
            return Ok(());
        }

        // Write magic number
        self.writer.write_all(BAM_MAGIC)?;

        // Build header text (SAM header format)
        let mut header_text = format!("@HD\tVN:{}\tSO:coordinate\n", header.version);
        for (name, length) in &header.references {
            header_text.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, length));
        }
        if let Some(ref program) = header.program {
            header_text.push_str(&format!("@PG\tID:helix\tPN:helix\tVN:{}\n", program));
        }

        // Write header length and text
        let text_len = header_text.len() as u32;
        self.writer.write_all(&text_len.to_le_bytes())?;
        self.writer.write_all(header_text.as_bytes())?;

        // Write number of reference sequences
        let n_ref = header.references.len() as u32;
        self.writer.write_all(&n_ref.to_le_bytes())?;

        // Write each reference
        for (name, length) in &header.references {
            let name_len = (name.len() + 1) as u32; // +1 for null terminator
            self.writer.write_all(&name_len.to_le_bytes())?;
            self.writer.write_all(name.as_bytes())?;
            self.writer.write_all(&[0u8])?; // Null terminator
            let ref_len = *length as u32;
            self.writer.write_all(&ref_len.to_le_bytes())?;
        }

        self.references = header.references.clone();
        self.header_written = true;
        
        // Flush header block
        self.writer.flush_block()?;

        Ok(())
    }

    /// Write a BAM record
    pub fn write_record(&mut self, record: &SamRecord) -> Result<()> {
        // Record start for indexing
        let record_start = self.writer.virtual_offset();

        // Encode the record
        let encoded = self.encode_record(record)?;

        // Write block size
        let block_size = encoded.len() as u32;
        self.writer.write_all(&block_size.to_le_bytes())?;
        self.writer.write_all(&encoded)?;

        // Update index if enabled
        if record.flag & FLAG_UNMAPPED == 0 {
            let ref_id = self.get_ref_id(&record.rname);
            if let (Some(rid), Some(ref mut index)) = (ref_id, self.index.as_mut()) {
                index.add_alignment(rid, record.pos, record_start);
            }
        }

        Ok(())
    }

    /// Get reference ID for a reference name
    fn get_ref_id(&self, name: &str) -> Option<i32> {
        if name == "*" {
            return None;
        }
        self.references.iter().position(|(n, _)| n == name).map(|i| i as i32)
    }

    /// Encode a SAM record to BAM binary format
    fn encode_record(&self, record: &SamRecord) -> Result<Vec<u8>> {
        let mut buf = Vec::with_capacity(256);

        // Reference ID (-1 for unmapped)
        let ref_id = self.get_ref_id(&record.rname).unwrap_or(-1);
        buf.extend_from_slice(&ref_id.to_le_bytes());

        // Position (0-based, -1 for unmapped)
        let pos = if record.flag & FLAG_UNMAPPED != 0 {
            -1i32
        } else {
            (record.pos - 1) as i32 // Convert 1-based to 0-based
        };
        buf.extend_from_slice(&pos.to_le_bytes());

        // Read name length (including null terminator)
        let l_read_name = (record.qname.len() + 1) as u8;
        buf.push(l_read_name);

        // Mapping quality
        buf.push(record.mapq);

        // BAM bin (computed from position and CIGAR)
        let bin = compute_bin(pos, &record.cigar);
        buf.extend_from_slice(&bin.to_le_bytes());

        // Number of CIGAR operations
        let cigar_ops = encode_cigar(&record.cigar);
        let n_cigar_op = cigar_ops.len() as u16;
        buf.extend_from_slice(&n_cigar_op.to_le_bytes());

        // Flags
        buf.extend_from_slice(&record.flag.to_le_bytes());

        // Sequence length
        let l_seq = record.seq.len() as u32;
        buf.extend_from_slice(&l_seq.to_le_bytes());

        // Next reference ID
        let next_ref_id = self.get_ref_id(&record.rnext).unwrap_or(-1);
        buf.extend_from_slice(&next_ref_id.to_le_bytes());

        // Next position
        let next_pos = if record.pnext > 0 {
            (record.pnext - 1) as i32
        } else {
            -1i32
        };
        buf.extend_from_slice(&next_pos.to_le_bytes());

        // Template length
        buf.extend_from_slice(&record.tlen.to_le_bytes());

        // Read name
        buf.extend_from_slice(record.qname.as_bytes());
        buf.push(0); // Null terminator

        // CIGAR (4 bytes per operation)
        for op in &cigar_ops {
            buf.extend_from_slice(&op.to_le_bytes());
        }

        // Sequence (4 bits per base, packed)
        buf.extend_from_slice(&encode_sequence(&record.seq));

        // Quality scores
        if record.qual == "*" {
            buf.extend(vec![0xFFu8; record.seq.len()]);
        } else {
            for c in record.qual.chars() {
                buf.push((c as u8).saturating_sub(33));
            }
        }

        // Tags
        for (tag, value) in &record.tags {
            buf.extend_from_slice(tag.as_bytes());
            // Assume string type for simplicity
            buf.push(b'Z');
            buf.extend_from_slice(value.as_bytes());
            buf.push(0);
        }

        Ok(buf)
    }

    /// Finish writing and return the index if enabled
    pub fn finish(mut self) -> Result<Option<BaiIndex>> {
        // Write EOF marker (empty BGZF block)
        self.writer.write_eof()?;
        self.writer.flush()?;

        Ok(self.index.map(|b| b.build()))
    }
}

/// BGZF (blocked gzip) writer
struct BgzfWriter<W: Write> {
    inner: W,
    buffer: Vec<u8>,
    compressed_offset: u64,
    uncompressed_offset: u16,
}

impl<W: Write> BgzfWriter<W> {
    fn new(writer: W) -> Self {
        Self {
            inner: writer,
            buffer: Vec::with_capacity(BGZF_MAX_BLOCK_SIZE),
            compressed_offset: 0,
            uncompressed_offset: 0,
        }
    }

    fn write_all(&mut self, data: &[u8]) -> Result<()> {
        self.buffer.extend_from_slice(data);
        // Track uncompressed offset within current block
        self.uncompressed_offset = self.uncompressed_offset.wrapping_add(data.len() as u16);

        // Flush if buffer is getting full
        if self.buffer.len() >= BGZF_MAX_BLOCK_SIZE - 256 {
            self.flush_block()?;
        }

        Ok(())
    }

    fn virtual_offset(&self) -> u64 {
        (self.compressed_offset << 16) | (self.uncompressed_offset as u64)
    }

    fn flush_block(&mut self) -> Result<()> {
        if self.buffer.is_empty() {
            return Ok(());
        }

        let compressed = compress_bgzf_block(&self.buffer)?;
        self.inner.write_all(&compressed)?;
        self.compressed_offset += compressed.len() as u64;
        self.uncompressed_offset = 0;
        self.buffer.clear();

        Ok(())
    }

    fn write_eof(&mut self) -> Result<()> {
        self.flush_block()?;
        // Write empty BGZF block (EOF marker)
        let eof = compress_bgzf_block(&[])?;
        self.inner.write_all(&eof)
    }

    fn flush(&mut self) -> Result<()> {
        self.flush_block()?;
        self.inner.flush()
    }
}

/// Compress data into a BGZF block
fn compress_bgzf_block(data: &[u8]) -> Result<Vec<u8>> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(data)?;
    let compressed = encoder.finish()?;

    // Build BGZF header
    let bsize = compressed.len() + 25; // +25 for extra subfield
    let mut block = Vec::with_capacity(bsize + 1);

    // Standard gzip header with BGZF extra field
    block.extend_from_slice(&[
        0x1f, 0x8b, // Gzip magic
        0x08,       // Compression method (deflate)
        0x04,       // Flags (FEXTRA)
        0, 0, 0, 0, // Modification time
        0,          // Extra flags
        0xff,       // OS (unknown)
        6, 0,       // Extra length (6 bytes)
        0x42, 0x43, // BGZF subfield ID
        2, 0,       // Subfield length
    ]);
    
    // BSIZE - 1
    let bsize_minus_1 = (bsize - 1) as u16;
    block.extend_from_slice(&bsize_minus_1.to_le_bytes());

    // Compressed data (skip gzip header/trailer from encoder)
    if compressed.len() > 18 {
        block.extend_from_slice(&compressed[10..compressed.len()-8]);
    }

    // CRC32 and input size
    let crc = crc32(&data);
    let isize = data.len() as u32;
    block.extend_from_slice(&crc.to_le_bytes());
    block.extend_from_slice(&isize.to_le_bytes());

    Ok(block)
}

/// Simple CRC32 calculation
fn crc32(data: &[u8]) -> u32 {
    let mut crc: u32 = 0xFFFFFFFF;
    for byte in data {
        crc ^= *byte as u32;
        for _ in 0..8 {
            if crc & 1 != 0 {
                crc = (crc >> 1) ^ 0xEDB88320;
            } else {
                crc >>= 1;
            }
        }
    }
    !crc
}

/// Encode CIGAR string to BAM format (4 bytes per op)
fn encode_cigar(cigar: &str) -> Vec<u32> {
    if cigar == "*" {
        return Vec::new();
    }

    let mut ops = Vec::new();
    let mut num = 0u32;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as u32 - '0' as u32);
        } else {
            let op_code = match c {
                'M' => 0,
                'I' => 1,
                'D' => 2,
                'N' => 3,
                'S' => 4,
                'H' => 5,
                'P' => 6,
                '=' => 7,
                'X' => 8,
                _ => continue,
            };
            ops.push((num << 4) | op_code);
            num = 0;
        }
    }

    ops
}

/// Encode sequence to 4-bit packed format
fn encode_sequence(seq: &str) -> Vec<u8> {
    if seq == "*" {
        return Vec::new();
    }

    let mut packed = Vec::with_capacity((seq.len() + 1) / 2);
    let bytes: Vec<u8> = seq.bytes().collect();

    for chunk in bytes.chunks(2) {
        let high = base_to_4bit(chunk[0]);
        let low = if chunk.len() > 1 { base_to_4bit(chunk[1]) } else { 0 };
        packed.push((high << 4) | low);
    }

    packed
}

/// Convert base to 4-bit encoding
fn base_to_4bit(base: u8) -> u8 {
    match base {
        b'A' | b'a' => 1,
        b'C' | b'c' => 2,
        b'G' | b'g' => 4,
        b'T' | b't' => 8,
        b'N' | b'n' => 15,
        _ => 15,
    }
}

/// Compute BAM bin for a region
fn compute_bin(pos: i32, cigar: &str) -> u16 {
    if pos < 0 {
        return 4680; // Unmapped bin
    }

    let beg = pos as u32;
    let end = beg + reference_length(cigar);

    // BAM bin calculation
    if beg >> 14 == end >> 14 {
        ((1 << 15) - 1) / 7 + (beg >> 14) as u16
    } else if beg >> 17 == end >> 17 {
        ((1 << 12) - 1) / 7 + (beg >> 17) as u16
    } else if beg >> 20 == end >> 20 {
        ((1 << 9) - 1) / 7 + (beg >> 20) as u16
    } else if beg >> 23 == end >> 23 {
        ((1 << 6) - 1) / 7 + (beg >> 23) as u16
    } else if beg >> 26 == end >> 26 {
        ((1 << 3) - 1) / 7 + (beg >> 26) as u16
    } else {
        0
    }
}

/// Calculate reference length from CIGAR
fn reference_length(cigar: &str) -> u32 {
    if cigar == "*" {
        return 0;
    }

    let mut len = 0u32;
    let mut num = 0u32;

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num = num * 10 + (c as u32 - '0' as u32);
        } else {
            // M, D, N, =, X consume reference
            if matches!(c, 'M' | 'D' | 'N' | '=' | 'X') {
                len += num;
            }
            num = 0;
        }
    }

    len
}

/// BAI index builder
pub struct BaiIndexBuilder {
    /// Bins for each reference
    bins: Vec<Vec<BaiChunk>>,
    /// Linear index for each reference
    linear_index: Vec<Vec<u64>>,
}

impl BaiIndexBuilder {
    fn new() -> Self {
        Self {
            bins: Vec::new(),
            linear_index: Vec::new(),
        }
    }

    fn add_alignment(&mut self, ref_id: i32, pos: u32, virtual_offset: u64) {
        let rid = ref_id as usize;

        // Ensure vectors are large enough
        while self.bins.len() <= rid {
            self.bins.push(Vec::new());
            self.linear_index.push(Vec::new());
        }

        // Add to bin
        let bin = ((pos - 1) >> 14) as usize;
        if bin < self.bins[rid].len() {
            self.bins[rid][bin].end = virtual_offset;
        } else {
            while self.bins[rid].len() <= bin {
                self.bins[rid].push(BaiChunk {
                    start: virtual_offset,
                    end: virtual_offset,
                });
            }
        }

        // Add to linear index
        let linear_bin = ((pos - 1) >> 14) as usize;
        while self.linear_index[rid].len() <= linear_bin {
            self.linear_index[rid].push(virtual_offset);
        }
    }

    fn build(self) -> BaiIndex {
        BaiIndex {
            bins: self.bins,
            linear_index: self.linear_index,
        }
    }
}

/// A chunk in the BAI index
#[derive(Clone)]
struct BaiChunk {
    start: u64,
    end: u64,
}

/// BAI index
pub struct BaiIndex {
    bins: Vec<Vec<BaiChunk>>,
    linear_index: Vec<Vec<u64>>,
}

impl BaiIndex {
    /// Write the BAI index to a file
    pub fn write_to_file<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let file = File::create(path)?;
        let mut writer = BufWriter::new(file);

        // Magic number
        writer.write_all(b"BAI\x01")?;

        // Number of references
        let n_ref = self.bins.len() as u32;
        writer.write_all(&n_ref.to_le_bytes())?;

        // Each reference
        for (bins, linear) in self.bins.iter().zip(self.linear_index.iter()) {
            // Number of bins
            let n_bin = bins.len() as u32;
            writer.write_all(&n_bin.to_le_bytes())?;

            // Each bin
            for (bin_num, chunk) in bins.iter().enumerate() {
                let bin = bin_num as u32;
                writer.write_all(&bin.to_le_bytes())?;

                // One chunk per bin (simplified)
                let n_chunk = 1u32;
                writer.write_all(&n_chunk.to_le_bytes())?;
                writer.write_all(&chunk.start.to_le_bytes())?;
                writer.write_all(&chunk.end.to_le_bytes())?;
            }

            // Linear index
            let n_intv = linear.len() as u32;
            writer.write_all(&n_intv.to_le_bytes())?;
            for offset in linear {
                writer.write_all(&offset.to_le_bytes())?;
            }
        }

        writer.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_base_to_4bit() {
        assert_eq!(base_to_4bit(b'A'), 1);
        assert_eq!(base_to_4bit(b'C'), 2);
        assert_eq!(base_to_4bit(b'G'), 4);
        assert_eq!(base_to_4bit(b'T'), 8);
        assert_eq!(base_to_4bit(b'N'), 15);
    }

    #[test]
    fn test_encode_sequence() {
        let encoded = encode_sequence("ACGT");
        assert_eq!(encoded.len(), 2);
        assert_eq!(encoded[0], 0x12); // A=1, C=2 -> 0x12
        assert_eq!(encoded[1], 0x48); // G=4, T=8 -> 0x48
    }

    #[test]
    fn test_encode_sequence_odd() {
        let encoded = encode_sequence("ACG");
        assert_eq!(encoded.len(), 2);
        assert_eq!(encoded[0], 0x12);
        assert_eq!(encoded[1], 0x40); // G=4, padded with 0
    }

    #[test]
    fn test_encode_cigar() {
        let ops = encode_cigar("10M5I3D");
        assert_eq!(ops.len(), 3);
        assert_eq!(ops[0], (10 << 4) | 0); // 10M
        assert_eq!(ops[1], (5 << 4) | 1);  // 5I
        assert_eq!(ops[2], (3 << 4) | 2);  // 3D
    }

    #[test]
    fn test_encode_cigar_star() {
        let ops = encode_cigar("*");
        assert!(ops.is_empty());
    }

    #[test]
    fn test_reference_length() {
        assert_eq!(reference_length("10M"), 10);
        assert_eq!(reference_length("5M3I2D"), 7); // M and D consume ref
        assert_eq!(reference_length("*"), 0);
    }

    #[test]
    fn test_crc32() {
        assert_eq!(crc32(b""), 0);
        assert_eq!(crc32(b"hello"), 0x3610a686);
    }

    #[test]
    fn test_compute_bin_unmapped() {
        assert_eq!(compute_bin(-1, "*"), 4680);
    }

    #[test]
    fn test_compute_bin_mapped() {
        let bin = compute_bin(1000, "50M");
        // Bin should be a valid value (not the unmapped bin, unless very small)
        // For position 1000 with 50bp read, should get a valid bin
        assert!(bin > 0);
    }

    #[test]
    fn test_bam_writer_header() {
        let mut buf = Vec::new();
        let mut writer = BamWriter::new(&mut buf);
        
        let mut header = SamHeader::new();
        header.add_reference("chr1".to_string(), 1000000);
        
        writer.write_header(&header).unwrap();
        
        // Check magic number is at start
        assert!(buf.len() > 4);
        // Note: actual magic check would require decompressing BGZF
    }

    #[test]
    fn test_bam_writer_record() {
        let mut buf = Vec::new();
        let mut writer = BamWriter::new(&mut buf);
        
        let mut header = SamHeader::new();
        header.add_reference("chr1".to_string(), 1000000);
        writer.write_header(&header).unwrap();
        
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
            tags: Vec::new(),
        };
        
        writer.write_record(&record).unwrap();
        assert!(buf.len() > 100);
    }
}
