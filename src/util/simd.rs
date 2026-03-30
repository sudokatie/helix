/// SIMD capability levels
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SimdLevel {
    /// No SIMD acceleration available
    None,
    /// x86-64 AVX2 (256-bit vectors)
    Avx2,
    /// ARM64 NEON (128-bit vectors)
    Neon,
}

impl SimdLevel {
    /// Human-readable name
    pub fn name(&self) -> &'static str {
        match self {
            SimdLevel::None => "scalar",
            SimdLevel::Avx2 => "AVX2",
            SimdLevel::Neon => "NEON",
        }
    }
}

/// Detect the best available SIMD level at runtime
pub fn detect_simd() -> SimdLevel {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("avx2") {
            return SimdLevel::Avx2;
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // NEON is always available on aarch64
        return SimdLevel::Neon;
    }

    SimdLevel::None
}

/// Check if AVX2 is available (x86-64 only)
#[cfg(target_arch = "x86_64")]
pub fn is_avx2_available() -> bool {
    is_x86_feature_detected!("avx2")
}

#[cfg(not(target_arch = "x86_64"))]
pub fn is_avx2_available() -> bool {
    false
}

/// Check if NEON is available (always true on aarch64)
#[cfg(target_arch = "aarch64")]
pub fn is_neon_available() -> bool {
    true
}

#[cfg(not(target_arch = "aarch64"))]
pub fn is_neon_available() -> bool {
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_simd() {
        let level = detect_simd();
        // Should detect something on modern hardware
        println!("Detected SIMD level: {:?} ({})", level, level.name());
        
        // Just verify it returns a valid level
        assert!(matches!(
            level,
            SimdLevel::None | SimdLevel::Avx2 | SimdLevel::Neon
        ));
    }

    #[test]
    fn test_simd_level_name() {
        assert_eq!(SimdLevel::None.name(), "scalar");
        assert_eq!(SimdLevel::Avx2.name(), "AVX2");
        assert_eq!(SimdLevel::Neon.name(), "NEON");
    }

    #[test]
    fn test_platform_detection() {
        #[cfg(target_arch = "x86_64")]
        {
            // On x86-64, AVX2 detection should work
            let _ = is_avx2_available();
        }

        #[cfg(target_arch = "aarch64")]
        {
            // On aarch64, NEON should always be available
            assert!(is_neon_available());
        }
    }
}
