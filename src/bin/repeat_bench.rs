use gaussian_crypto::{field_mul, rand_nonzero_below_p, U256};
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let runs = 8usize;
    let iters = 20_000u64;
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
        println!("run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;

    println!("summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}
