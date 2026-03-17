use gaussian_crypto::{field_inv, field_mul, P, U256};
use std::hint::black_box;
use std::time::Instant;

#[derive(Clone, Copy)]
enum InvVariant {
    PowBaseline,
    BatchAll,
    BatchChunks32,
}

impl InvVariant {
    fn name(&self) -> &'static str {
        match self {
            InvVariant::PowBaseline => "pow_baseline",
            InvVariant::BatchAll => "batch_all",
            InvVariant::BatchChunks32 => "batch_chunks_32",
        }
    }
}

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *seed
}

fn rand_u256(seed: &mut u64) -> U256 {
    [lcg(seed), lcg(seed), lcg(seed), lcg(seed)]
}

fn cmp_u256(a: &U256, b: &U256) -> core::cmp::Ordering {
    for i in (0..4).rev() {
        if a[i] < b[i] {
            return core::cmp::Ordering::Less;
        } else if a[i] > b[i] {
            return core::cmp::Ordering::Greater;
        }
    }
    core::cmp::Ordering::Equal
}

fn sub_u256(a: U256, b: U256) -> (U256, u64) {
    let mut out = [0u64; 4];
    let mut borrow = 0u64;
    for i in 0..4 {
        let (r1, b1) = a[i].overflowing_sub(b[i]);
        let (r2, b2) = r1.overflowing_sub(borrow);
        out[i] = r2;
        borrow = (b1 as u64) | (b2 as u64);
    }
    (out, borrow)
}

fn reduce_once(x: U256) -> U256 {
    if cmp_u256(&x, &P) != core::cmp::Ordering::Less {
        let (y, _) = sub_u256(x, P);
        y
    } else {
        x
    }
}

fn rand_nonzero_below_p(seed: &mut u64) -> U256 {
    loop {
        let x = reduce_once(rand_u256(seed));
        if x != [0, 0, 0, 0] {
            return x;
        }
    }
}

fn inv_pow_baseline(a: U256) -> U256 {
    field_inv(a)
}

fn batch_invert(values: &[U256]) -> Vec<U256> {
    let n = values.len();
    if n == 0 {
        return vec![];
    }

    let mut prefix = vec![[0u64; 4]; n];
    let mut acc = [1u64, 0, 0, 0];

    for (i, &v) in values.iter().enumerate() {
        prefix[i] = acc;
        acc = field_mul(acc, v);
    }

    let mut inv_acc = field_inv(acc);
    let mut out = vec![[0u64; 4]; n];

    for i in (0..n).rev() {
        out[i] = field_mul(inv_acc, prefix[i]);
        inv_acc = field_mul(inv_acc, values[i]);
    }

    out
}

fn invert_variant(variant: InvVariant, inputs: &[U256]) -> Vec<U256> {
    match variant {
        InvVariant::PowBaseline => inputs.iter().copied().map(inv_pow_baseline).collect(),
        InvVariant::BatchAll => batch_invert(inputs),
        InvVariant::BatchChunks32 => {
            let mut out = Vec::with_capacity(inputs.len());
            for chunk in inputs.chunks(32) {
                out.extend(batch_invert(chunk));
            }
            out
        }
    }
}

fn verify_inverses(inputs: &[U256], inverses: &[U256]) {
    assert_eq!(inputs.len(), inverses.len());
    for (&a, &inv) in inputs.iter().zip(inverses.iter()) {
        let one = field_mul(a, inv);
        assert_eq!(one, [1, 0, 0, 0], "inverse failed for {:?}", a);
    }
}

fn bench_mul(iters: u64) {
    let runs = 8usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xA11CE5EED1234567u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            let b = rand_nonzero_below_p(&mut seed);
            let c = field_mul(black_box(a), black_box(b));
            checksum ^= black_box(c[0]);
        }

        let dur = start.elapsed();
        println!("mul run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("mul summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_single_inv(variant: InvVariant, iters: usize) {
    let runs = 5usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xBEEFFEED12345678u64 ^ run as u64;
        let inputs: Vec<U256> = (0..iters).map(|_| rand_nonzero_below_p(&mut seed)).collect();

        let start = Instant::now();
        let inverses = invert_variant(variant, &inputs);
        let dur = start.elapsed();

        verify_inverses(&inputs, &inverses);
        let checksum = inverses.iter().fold(0u64, |acc, x| acc ^ x[0]);

        println!("{} single_inv run {}: {:?}, checksum={}", variant.name(), run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("{} single_inv summary_sec min={:.9} mean={:.9} max={:.9}", variant.name(), min, mean, max);
}

fn bench_batch_inv(variant: InvVariant, total_inputs: usize) {
    let runs = 5usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0x1234AAAABBBBCCCCu64 ^ run as u64;
        let inputs: Vec<U256> = (0..total_inputs).map(|_| rand_nonzero_below_p(&mut seed)).collect();

        let start = Instant::now();
        let inverses = invert_variant(variant, &inputs);
        let dur = start.elapsed();

        verify_inverses(&inputs, &inverses);
        let checksum = inverses.iter().fold(0u64, |acc, x| acc ^ x[0]);

        println!("{} batch_inv run {}: {:?}, checksum={}", variant.name(), run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("{} batch_inv summary_sec min={:.9} mean={:.9} max={:.9}", variant.name(), min, mean, max);
}

fn main() {
    let variants = [
        InvVariant::PowBaseline,
        InvVariant::BatchAll,
        InvVariant::BatchChunks32,
    ];

    println!("== tournament start ==");
    println!("spot check correctness first");

    let mut seed = 0xDEADBEEF12345678u64;
    let spot_inputs: Vec<U256> = (0..64).map(|_| rand_nonzero_below_p(&mut seed)).collect();
    for v in variants {
        let out = invert_variant(v, &spot_inputs);
        verify_inverses(&spot_inputs, &out);
        println!("{} correctness: PASS", v.name());
    }

    println!("\n== field mul benchmark ==");
    bench_mul(20_000);

    println!("\n== single inversion benchmark (2k inversions) ==");
    for v in variants {
        bench_single_inv(v, 2_000);
    }

    println!("\n== batch inversion benchmark (8k inversions total) ==");
    for v in variants {
        bench_batch_inv(v, 8_000);
    }

    println!("\n== tournament done ==");
}
