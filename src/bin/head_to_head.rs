use gaussian_crypto::{field_mul, rand_nonzero_below_p, U256};
use k256::Scalar;
use std::hint::black_box;
use std::time::Instant;

fn bench_mine(iters: u64) {
    let runs = 8usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xCAFEBABE12345678u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a: U256 = rand_nonzero_below_p(&mut seed);
            let b: U256 = rand_nonzero_below_p(&mut seed);
            let c = field_mul(black_box(a), black_box(b));
            checksum ^= black_box(c[0]);
        }

        let dur = start.elapsed();
        println!("mine run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("mine summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_k256(iters: u64) {
    let runs = 8usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let start = Instant::now();
        let mut checksum = 0u64;

        for i in 0..iters {
            let a = Scalar::from(i.wrapping_add((run as u64) << 32));
            let b = Scalar::from(i.wrapping_mul(0x9E3779B185EBCA87).wrapping_add(17));
            let c = black_box(a) * black_box(b);

            let bytes = c.to_bytes();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[24..32]);
            checksum ^= u64::from_be_bytes(low);
        }

        let dur = start.elapsed();
        println!("k256 run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("k256 summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn main() {
    let iters = 20_000u64;
    println!("== head-to-head benchmark: {} iterations per run ==", iters);
    bench_mine(iters);
    bench_k256(iters);
}
