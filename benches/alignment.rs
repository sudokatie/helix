/// Benchmarks for helix alignment

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};

use helix::align::{align, smith_waterman_scalar, AlignmentConfig};
use helix::seq::encode_sequence;

/// Generate a random-ish sequence of given length
fn generate_sequence(len: usize, seed: u64) -> Vec<u8> {
    let bases = [0u8, 1, 2, 3]; // A, C, G, T encoded
    let mut seq = Vec::with_capacity(len);
    let mut state = seed;

    for _ in 0..len {
        // Simple LCG
        state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        seq.push(bases[(state >> 32) as usize % 4]);
    }

    seq
}

fn bench_smith_waterman_scalar(c: &mut Criterion) {
    let config = AlignmentConfig::default();

    let mut group = c.benchmark_group("sw_scalar");

    for size in [50, 100, 200, 500].iter() {
        let query = generate_sequence(*size, 12345);
        let target = generate_sequence(*size, 67890);

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| smith_waterman_scalar(black_box(&query), black_box(&target), &config))
        });
    }

    group.finish();
}

fn bench_smith_waterman_dispatch(c: &mut Criterion) {
    let config = AlignmentConfig::default();

    let mut group = c.benchmark_group("sw_dispatch");

    for size in [50, 100, 200, 500].iter() {
        let query = generate_sequence(*size, 12345);
        let target = generate_sequence(*size, 67890);

        group.bench_with_input(BenchmarkId::from_parameter(size), size, |b, _| {
            b.iter(|| align(black_box(&query), black_box(&target), &config))
        });
    }

    group.finish();
}

fn bench_exact_match(c: &mut Criterion) {
    let config = AlignmentConfig::default();

    // Exact match (best case)
    let seq = generate_sequence(200, 11111);

    c.bench_function("sw_exact_200bp", |b| {
        b.iter(|| smith_waterman_scalar(black_box(&seq), black_box(&seq), &config))
    });
}

fn bench_with_errors(c: &mut Criterion) {
    let config = AlignmentConfig::default();

    // ~5% error rate
    let query = generate_sequence(200, 22222);
    let mut target = query.clone();
    for i in (0..target.len()).step_by(20) {
        target[i] = (target[i] + 1) % 4; // Introduce mismatch
    }

    c.bench_function("sw_5pct_error_200bp", |b| {
        b.iter(|| smith_waterman_scalar(black_box(&query), black_box(&target), &config))
    });
}

criterion_group!(
    benches,
    bench_smith_waterman_scalar,
    bench_smith_waterman_dispatch,
    bench_exact_match,
    bench_with_errors
);
criterion_main!(benches);
