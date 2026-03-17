use gaussian_crypto::{field_inv, field_mul, field_square, P, U256};
use k256::{FieldElement, ProjectivePoint};
use k256::elliptic_curve::BatchNormalize;
use k256::elliptic_curve::sec1::ToEncodedPoint;
use std::hint::black_box;
use std::time::Instant;

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

fn u256_to_be32(x: U256) -> [u8; 32] {
    let mut out = [0u8; 32];
    for i in 0..4 {
        out[(3 - i) * 8..(4 - i) * 8].copy_from_slice(&x[i].to_be_bytes());
    }
    out
}

#[allow(dead_code)]
fn to_fe(x: U256) -> FieldElement {
    let bytes = u256_to_be32(x);
    let ct = FieldElement::from_bytes((&bytes).into());
    Option::<FieldElement>::from(ct).expect("valid field element")
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

fn bench_mine(points: usize) {
    let runs = 6usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut seed = 0xABCD1234EF998877u64 ^ run as u64;
        let xs: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();
        let ys: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();
        let zs: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();

        let start = Instant::now();

        let zinv = batch_invert(&zs);
        let mut checksum = 0u64;

        for i in 0..points {
            let z2 = field_square(zinv[i]);
            let z3 = field_mul(z2, zinv[i]);
            let ax = field_mul(xs[i], z2);
            let ay = field_mul(ys[i], z3);
            checksum ^= black_box(ax[0] ^ ay[0]);
        }

        let dur = start.elapsed();
        println!("mine_norm run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("mine_norm summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_k256(points: usize) {
    let runs = 6usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut pts = Vec::with_capacity(points);

        for i in 0..points {
            let s = (i as u64).wrapping_add(1).wrapping_add((run as u64) << 32);
            let p = ProjectivePoint::GENERATOR * black_box(k256::Scalar::from(s));
            pts.push(p);
        }

        let start = Instant::now();
        let out = ProjectivePoint::batch_normalize(pts.as_slice());
        let dur = start.elapsed();

        let mut checksum = 0u64;
        for p in &out {
            let enc = p.to_encoded_point(false);
            let bytes = enc.as_bytes();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[bytes.len() - 8..]);
            checksum ^= u64::from_be_bytes(low);
        }

        println!("k256_norm run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("k256_norm summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn main() {
    let points = 4096usize;
    println!("== real workload: batch normalization-style benchmark ==");
    println!("points per run: {}", points);
    bench_mine(points);
    bench_k256(points);
}
