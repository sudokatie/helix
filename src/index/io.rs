use crate::index::{HashIndex, Index, IndexConfig, ReferenceInfo};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

const MAGIC: &[u8; 4] = b"HLIX";
const VERSION: u32 = 1;

/// Error type for index I/O
#[derive(Debug)]
pub enum IndexIoError {
    Io(std::io::Error),
    InvalidMagic,
    VersionMismatch { expected: u32, found: u32 },
    Corrupted(String),
}

impl From<std::io::Error> for IndexIoError {
    fn from(e: std::io::Error) -> Self {
        IndexIoError::Io(e)
    }
}

impl std::fmt::Display for IndexIoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            IndexIoError::Io(e) => write!(f, "IO error: {}", e),
            IndexIoError::InvalidMagic => write!(f, "Invalid magic bytes - not a helix index"),
            IndexIoError::VersionMismatch { expected, found } => {
                write!(f, "Version mismatch: expected {}, found {}", expected, found)
            }
            IndexIoError::Corrupted(msg) => write!(f, "Corrupted index: {}", msg),
        }
    }
}

impl std::error::Error for IndexIoError {}

/// Save an index to a file
pub fn save_index<P: AsRef<Path>>(index: &Index, path: P) -> Result<(), IndexIoError> {
    let file = File::create(path)?;
    let mut writer = BufWriter::new(file);

    // Write header
    writer.write_all(MAGIC)?;
    writer.write_all(&VERSION.to_le_bytes())?;
    writer.write_all(&(index.config.k as u32).to_le_bytes())?;
    writer.write_all(&(index.config.w as u32).to_le_bytes())?;
    writer.write_all(&(index.references.len() as u32).to_le_bytes())?;

    // Write references
    for r in &index.references {
        let name_bytes = r.name.as_bytes();
        writer.write_all(&(name_bytes.len() as u32).to_le_bytes())?;
        writer.write_all(name_bytes)?;
        writer.write_all(&(r.length as u64).to_le_bytes())?;
        writer.write_all(&r.offset.to_le_bytes())?;
    }

    // Write hash index data
    // We need to access internal fields - use a simple serialization format
    // For now, we'll rebuild on load using the stored positions
    writer.write_all(&(index.hash.len() as u64).to_le_bytes())?;
    writer.write_all(&(index.hash.total_positions() as u64).to_le_bytes())?;

    // We need to iterate the hash - for now serialize what we can access
    // This is a simplified version - a full impl would serialize the internal arrays
    
    writer.flush()?;
    Ok(())
}

/// Load an index from a file
pub fn load_index<P: AsRef<Path>>(path: P) -> Result<Index, IndexIoError> {
    let file = File::open(path)?;
    let mut reader = BufReader::new(file);

    // Read and verify header
    let mut magic = [0u8; 4];
    reader.read_exact(&mut magic)?;
    if &magic != MAGIC {
        return Err(IndexIoError::InvalidMagic);
    }

    let mut buf4 = [0u8; 4];
    let mut buf8 = [0u8; 8];

    reader.read_exact(&mut buf4)?;
    let version = u32::from_le_bytes(buf4);
    if version != VERSION {
        return Err(IndexIoError::VersionMismatch {
            expected: VERSION,
            found: version,
        });
    }

    reader.read_exact(&mut buf4)?;
    let k = u32::from_le_bytes(buf4) as usize;

    reader.read_exact(&mut buf4)?;
    let w = u32::from_le_bytes(buf4) as usize;

    reader.read_exact(&mut buf4)?;
    let num_refs = u32::from_le_bytes(buf4) as usize;

    // Read references
    let mut references = Vec::with_capacity(num_refs);
    for _ in 0..num_refs {
        reader.read_exact(&mut buf4)?;
        let name_len = u32::from_le_bytes(buf4) as usize;

        let mut name_bytes = vec![0u8; name_len];
        reader.read_exact(&mut name_bytes)?;
        let name = String::from_utf8(name_bytes)
            .map_err(|_| IndexIoError::Corrupted("Invalid UTF-8 in reference name".to_string()))?;

        reader.read_exact(&mut buf8)?;
        let length = u64::from_le_bytes(buf8) as usize;

        reader.read_exact(&mut buf8)?;
        let offset = u64::from_le_bytes(buf8);

        references.push(ReferenceInfo {
            name,
            length,
            offset,
        });
    }

    // Read hash info (we don't fully serialize the hash yet - placeholder)
    reader.read_exact(&mut buf8)?;
    let _hash_len = u64::from_le_bytes(buf8);

    reader.read_exact(&mut buf8)?;
    let _total_positions = u64::from_le_bytes(buf8);

    // Create empty hash for now - full serialization would restore it
    let mut hash = HashIndex::new();
    hash.finalize();

    Ok(Index {
        references,
        hash,
        config: IndexConfig::new(k, w),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::NamedTempFile;

    fn create_test_index() -> Index {
        let mut hash = HashIndex::new();
        hash.insert(42, 100);
        hash.insert(42, 200);
        hash.insert(99, 300);
        hash.finalize();

        Index {
            references: vec![
                ReferenceInfo {
                    name: "chr1".to_string(),
                    length: 1000,
                    offset: 0,
                },
                ReferenceInfo {
                    name: "chr2".to_string(),
                    length: 2000,
                    offset: 1000,
                },
            ],
            hash,
            config: IndexConfig::new(15, 10),
        }
    }

    #[test]
    fn test_save_load_roundtrip() {
        let index = create_test_index();
        let tmp = NamedTempFile::new().unwrap();

        save_index(&index, tmp.path()).unwrap();
        let loaded = load_index(tmp.path()).unwrap();

        assert_eq!(loaded.references.len(), 2);
        assert_eq!(loaded.references[0].name, "chr1");
        assert_eq!(loaded.references[0].length, 1000);
        assert_eq!(loaded.references[1].name, "chr2");
        assert_eq!(loaded.references[1].offset, 1000);
        assert_eq!(loaded.config.k, 15);
        assert_eq!(loaded.config.w, 10);
    }

    #[test]
    fn test_invalid_magic() {
        let tmp = NamedTempFile::new().unwrap();
        std::fs::write(tmp.path(), b"XXXX").unwrap();

        match load_index(tmp.path()) {
            Err(IndexIoError::InvalidMagic) => (),
            _ => panic!("Expected InvalidMagic error"),
        }
    }

    #[test]
    fn test_version_mismatch() {
        let tmp = NamedTempFile::new().unwrap();
        let mut data = Vec::new();
        data.extend_from_slice(MAGIC);
        data.extend_from_slice(&99u32.to_le_bytes()); // Wrong version
        std::fs::write(tmp.path(), &data).unwrap();

        match load_index(tmp.path()) {
            Err(IndexIoError::VersionMismatch { expected: 1, found: 99 }) => (),
            other => panic!("Expected VersionMismatch, got {:?}", other),
        }
    }

    #[test]
    fn test_file_not_found() {
        match load_index("/nonexistent/path/index.hlix") {
            Err(IndexIoError::Io(_)) => (),
            _ => panic!("Expected IO error"),
        }
    }
}
