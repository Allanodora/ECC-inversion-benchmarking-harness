use gaussian_crypto::{field_inv, field_mul, field_square, P, U256};
use k256::{FieldElement, ProjectivePoint};
use k256::elliptic_curve::BatchNormalize;
use k256::elliptic_curve::sec1::ToEncodedPoint;
use std::hint::black_box;
use std::time::Instant;

#[derive(Clone, Copy)]
enum NormVariant {
    SingleInv,
    BatchAll,
    BatchChunks32,
    BatchChunks64,
}

impl NormVariant {
    fn name(&self) -> &'static str {
        match self {
            NormVariant::SingleInv => "single_inv",
            NormVariant::BatchAll => "batch_all",
            NormVariant::BatchChunks32 => "batch_chunks_32",
            NormVariant::BatchChunks64 => "batch_chunks_64",
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

fn u256_to_be32(x: U256) -> [u8; 32] {
    let mut out = [0u8; 32];
    for i in 0..4 {
        out[(3 - i) * 8..(4 - i) * 8].copy_from_slice(&x[i].to_be_bytes());
    }
    out
}

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

fn invert_many(variant: NormVariant, zs: &[U256]) -> Vec<U256> {
    match variant {
        NormVariant::SingleInv => zs.iter().copied().map(field_inv).collect(),
        NormVariant::BatchAll => batch_invert(zs),
        NormVariant::BatchChunks32 => {
            let mut out = Vec::with_capacity(zs.len());
            for chunk in zs.chunks(32) {
                out.extend(batch_invert(chunk));
            }
            out
        }
        NormVariant::BatchChunks64 => {
            let mut out = Vec::with_capacity(zs.len());
            for chunk in zs.chunks(64) {
                out.extend(batch_invert(chunk));
            }
            out
        }
    }
}

fn summarize(label: &str, values: &[f64]) {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let mean = sorted.iter().sum::<f64>() / sorted.len() as f64;
    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };
    println!(
        "{} min={:.9} median={:.9} mean={:.9} max={:.9}",
        label, min, median, mean, max
    );
}

fn bench_field_mul_mine(iters: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0x1111222233334444u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            let b = rand_nonzero_below_p(&mut seed);
            let c = field_mul(black_box(a), black_box(b));
            checksum ^= black_box(c[0]);
        }

        let dur = start.elapsed();
        println!("mine_mul run {:02}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize("mine_mul summary_sec", &results);
}

fn bench_field_mul_k256(iters: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0x1111222233334444u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = to_fe(rand_nonzero_below_p(&mut seed));
            let b = to_fe(rand_nonzero_below_p(&mut seed));
            let c = black_box(a) * black_box(b);
            let bytes = c.to_bytes();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[24..32]);
            checksum ^= u64::from_be_bytes(low);
        }

        let dur = start.elapsed();
        println!("k256_mul run {:02}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize("k256_mul summary_sec", &results);
}

fn bench_single_inv_mine(iters: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0x9999AAAABBBBCCCCu64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            let inv = field_inv(black_box(a));
            checksum ^= black_box(inv[0]);
        }

        let dur = start.elapsed();
        println!("mine_inv run {:02}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize("mine_inv summary_sec", &results);
}

fn bench_single_inv_k256(iters: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0x9999AAAABBBBCCCCu64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;

        for _ in 0..iters {
            let a = to_fe(rand_nonzero_below_p(&mut seed));
            let inv_ct = black_box(a).invert();
            let inv = Option::<FieldElement>::from(inv_ct).expect("nonzero invert");
            let bytes = inv.to_bytes();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[24..32]);
            checksum ^= u64::from_be_bytes(low);
        }

        let dur = start.elapsed();
        println!("k256_inv run {:02}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize("k256_inv summary_sec", &results);
}

fn bench_norm_mine(variant: NormVariant, points: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0xABCD1234EF998877u64 ^ run as u64;
        let xs: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();
        let ys: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();
        let zs: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();

        let start = Instant::now();
        let zinv = invert_many(variant, &zs);
        let mut checksum = 0u64;

        for i in 0..points {
            let z2 = field_square(zinv[i]);
            let z3 = field_mul(z2, zinv[i]);
            let ax = field_mul(xs[i], z2);
            let ay = field_mul(ys[i], z3);
            checksum ^= black_box(ax[0] ^ ay[0]);
        }

        let dur = start.elapsed();
        println!("{}_norm run {:02}: {:?}, checksum={}", variant.name(), run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize(&format!("{}_norm summary_sec", variant.name()), &results);
}

fn bench_norm_k256(points: usize, runs: usize) {
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

        println!("k256_norm run {:02}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize("k256_norm summary_sec", &results);
}

fn warmup() {
    let mut seed = 0x5555AAAACCCC9999u64;
    let mut sink = 0u64;
    for _ in 0..10_000 {
        let a = rand_nonzero_below_p(&mut seed);
        let b = rand_nonzero_below_p(&mut seed);
        let c = field_mul(a, b);
        let d = field_inv(a);
        sink ^= c[0] ^ d[0];
    }
    black_box(sink);
}

fn main() {
    let runs = 20usize;
    let mul_iters = 100_000usize;
    let inv_iters = 10_000usize;
    let norm_points = 8_192usize;

    println!("== production-style tournament ==");
    println!("warmup...");
    warmup();

    println!("\n== field multiplication ==");
    println!("iters per run: {}", mul_iters);
    bench_field_mul_mine(mul_iters, runs);
    bench_field_mul_k256(mul_iters, runs);

    println!("\n== single inversion ==");
    println!("iters per run: {}", inv_iters);
    bench_single_inv_mine(inv_iters, runs);
    bench_single_inv_k256(inv_iters, runs);

    println!("\n== real normalization workload ==");
    println!("points per run: {}", norm_points);
    bench_norm_mine(NormVariant::SingleInv, norm_points, runs);
    bench_norm_mine(NormVariant::BatchAll, norm_points, runs);
    bench_norm_mine(NormVariant::BatchChunks32, norm_points, runs);
    bench_norm_mine(NormVariant::BatchChunks64, norm_points, runs);
    bench_norm_k256(norm_points, runs);

    println!("\n== tournament complete ==");
}
