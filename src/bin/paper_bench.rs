use gaussian_crypto::{field_inv, field_mul, P, U256};
use k256::FieldElement;
use std::hint::black_box;
use std::time::Instant;

fn u256_to_be32(x: U256) -> [u8; 32] {
    let mut out = [0u8; 32];
    for i in 0..4 {
        out[(3 - i) * 8..(4 - i) * 8].copy_from_slice(&x[i].to_be_bytes());
    }
    out
}

fn be32_to_u256(x: [u8; 32]) -> U256 {
    let mut out = [0u64; 4];
    for i in 0..4 {
        let mut limb = [0u8; 8];
        limb.copy_from_slice(&x[(3 - i) * 8..(4 - i) * 8]);
        out[i] = u64::from_be_bytes(limb);
    }
    out
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

fn to_fe(x: U256) -> FieldElement {
    let bytes = u256_to_be32(x);
    let ct = FieldElement::from_bytes((&bytes).into());
    Option::<FieldElement>::from(ct).expect("valid field element")
}

fn from_fe(x: FieldElement) -> U256 {
    let bytes = x.to_bytes();
    let arr: [u8; 32] = bytes.into();
    be32_to_u256(arr)
}

fn bench_mul_mine(iters: u64) {
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
        println!("mine_mul run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("mine_mul summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_mul_k256(iters: u64) {
    let runs = 8usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xA11CE5EED1234567u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = to_fe(rand_nonzero_below_p(&mut seed));
            let b = to_fe(rand_nonzero_below_p(&mut seed));
            let c = black_box(a) * black_box(b);
            let u = from_fe(c);
            checksum ^= u[0];
        }

        let dur = start.elapsed();
        println!("k256_mul run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("k256_mul summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_inv_mine(iters: u64) {
    let runs = 5usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xBEEFFEED12345678u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            let inv = field_inv(black_box(a));
            checksum ^= black_box(inv[0]);
        }

        let dur = start.elapsed();
        println!("mine_inv run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("mine_inv summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_inv_k256(iters: u64) {
    let runs = 5usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xBEEFFEED12345678u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = to_fe(rand_nonzero_below_p(&mut seed));
            let inv_ct = black_box(a).invert();
            let inv = Option::<FieldElement>::from(inv_ct).expect("nonzero invert");
            let u = from_fe(inv);
            checksum ^= u[0];
        }

        let dur = start.elapsed();
        println!("k256_inv run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("k256_inv summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn main() {
    let mul_iters = 20_000u64;
    let inv_iters = 2_000u64;

    println!("== apples-to-apples field benchmark ==");
    println!("mul iters per run: {}", mul_iters);
    bench_mul_mine(mul_iters);
    bench_mul_k256(mul_iters);

    println!("inv iters per run: {}", inv_iters);
    bench_inv_mine(inv_iters);
    bench_inv_k256(inv_iters);
}
