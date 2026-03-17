use gaussian_crypto::{field_inv, field_mul, field_square, P, U256};
use k256::{FieldElement, ProjectivePoint};
use k256::elliptic_curve::BatchNormalize;
use k256::elliptic_curve::sec1::ToEncodedPoint;
use std::hint::black_box;
use std::time::Instant;

#[derive(Clone, Copy)]
enum NormVariant {
    SingleBaseline,
    SingleWindow4,
    BatchAllBaseline,
    BatchAllWindow4,
}

impl NormVariant {
    fn name(&self) -> &'static str {
        match self {
            NormVariant::SingleBaseline => "single_baseline",
            NormVariant::SingleWindow4 => "single_window4",
            NormVariant::BatchAllBaseline => "batch_all_baseline",
            NormVariant::BatchAllWindow4 => "batch_all_window4",
        }
    }
}

const ONE: U256 = [1, 0, 0, 0];
const P_MINUS_2: U256 = [
    0xFFFFFFFEFFFFFC2D,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
];

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

fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    }
}

fn summarize(label: &str, values: &[f64]) {
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min = sorted[0];
    let max = sorted[sorted.len() - 1];
    let mean = sorted.iter().sum::<f64>() / sorted.len() as f64;
    let med = median(values);
    println!(
        "{} min={:.9} median={:.9} mean={:.9} max={:.9}",
        label, min, med, mean, max
    );
}

fn inv_window4_v1(a: U256) -> U256 {
    if a == [0, 0, 0, 0] {
        return [0, 0, 0, 0];
    }

    let mut table = [[0u64; 4]; 16];
    table[0] = ONE;
    table[1] = a;
    for i in 2..16 {
        table[i] = field_mul(table[i - 1], a);
    }

    let exp = u256_to_be32(P_MINUS_2);
    let mut result = ONE;

    for byte in exp {
        let hi = (byte >> 4) as usize;
        let lo = (byte & 0x0f) as usize;

        for _ in 0..4 {
            result = field_square(result);
        }
        if hi != 0 {
            result = field_mul(result, table[hi]);
        }

        for _ in 0..4 {
            result = field_square(result);
        }
        if lo != 0 {
            result = field_mul(result, table[lo]);
        }
    }

    result
}

fn batch_invert_with<F>(values: &[U256], inv_fn: F) -> Vec<U256>
where
    F: Fn(U256) -> U256 + Copy,
{
    let n = values.len();
    if n == 0 {
        return vec![];
    }

    let mut prefix = vec![[0u64; 4]; n];
    let mut acc = ONE;

    for (i, &v) in values.iter().enumerate() {
        prefix[i] = acc;
        acc = field_mul(acc, v);
    }

    let mut inv_acc = inv_fn(acc);
    let mut out = vec![[0u64; 4]; n];

    for i in (0..n).rev() {
        out[i] = field_mul(inv_acc, prefix[i]);
        inv_acc = field_mul(inv_acc, values[i]);
    }

    out
}

fn invert_many(variant: NormVariant, zs: &[U256]) -> Vec<U256> {
    match variant {
        NormVariant::SingleBaseline => zs.iter().copied().map(field_inv).collect(),
        NormVariant::SingleWindow4 => zs.iter().copied().map(inv_window4_v1).collect(),
        NormVariant::BatchAllBaseline => batch_invert_with(zs, field_inv),
        NormVariant::BatchAllWindow4 => batch_invert_with(zs, inv_window4_v1),
    }
}

fn warmup() {
    let mut seed = 0x5555AAAACCCC9999u64;
    let mut sink = 0u64;
    for _ in 0..20_000 {
        let a = rand_nonzero_below_p(&mut seed);
        let b = rand_nonzero_below_p(&mut seed);
        let c = field_mul(a, b);
        let d = inv_window4_v1(a);
        sink ^= c[0] ^ d[0];
    }
    black_box(sink);
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

fn bench_single_inv_variant(label: &str, inv_fn: fn(U256) -> U256, iters: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0x9999AAAABBBBCCCCu64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;
        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            let inv = inv_fn(black_box(a));
            checksum ^= black_box(inv[0]);
        }
        let dur = start.elapsed();
        println!("{} run {:02}: {:?}, checksum={}", label, run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize(&format!("{} summary_sec", label), &results);
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

fn bench_norm_variant(variant: NormVariant, points: usize, runs: usize) {
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

fn scaling_sweep() {
    let sizes = [1024usize, 4096, 8192, 32768];
    let runs = 5usize;

    println!("\n== scaling sweep ==");
    for &points in &sizes {
        println!("points={}", points);
        bench_norm_variant(NormVariant::BatchAllBaseline, points, runs);
        bench_norm_variant(NormVariant::BatchAllWindow4, points, runs);
        bench_norm_k256(points, runs);
        println!();
    }
}

fn timing_sanity(inv_fn: fn(U256) -> U256, label: &str) {
    let runs = 10usize;
    let iters = 5000usize;
    let fixed = [0x1234567890ABCDEF, 0x0FEDCBA987654321, 1, 2];
    let mut fixed_times = Vec::with_capacity(runs);
    let mut rand_times = Vec::with_capacity(runs);

    for run in 0..runs {
        let start = Instant::now();
        let mut checksum = 0u64;
        for _ in 0..iters {
            let out = inv_fn(black_box(fixed));
            checksum ^= out[0];
        }
        let dur = start.elapsed();
        println!("{}_fixed run {:02}: {:?}, checksum={}", label, run + 1, dur, checksum);
        fixed_times.push(dur.as_secs_f64());
    }

    for run in 0..runs {
        let mut seed = 0xCAFED00D12345678u64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;
        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            let out = inv_fn(black_box(a));
            checksum ^= out[0];
        }
        let dur = start.elapsed();
        println!("{}_random run {:02}: {:?}, checksum={}", label, run + 1, dur, checksum);
        rand_times.push(dur.as_secs_f64());
    }

    summarize(&format!("{}_fixed summary_sec", label), &fixed_times);
    summarize(&format!("{}_random summary_sec", label), &rand_times);
}

fn correctness_spot() {
    let mut seed = 0xDEADBEEF12345678u64;
    for _ in 0..2000 {
        let a = rand_nonzero_below_p(&mut seed);
        let inv0 = field_inv(a);
        let inv1 = inv_window4_v1(a);
        let one0 = field_mul(a, inv0);
        let one1 = field_mul(a, inv1);
        assert_eq!(one0, ONE, "baseline inverse failed");
        assert_eq!(one1, ONE, "window4 inverse failed");
    }
    println!("correctness_spot: PASS");
}

fn main() {
    let runs = 20usize;
    let mul_iters = 100_000usize;
    let inv_iters = 10_000usize;
    let norm_points = 8_192usize;

    println!("== phase 2 compare ==");
    println!("warmup...");
    warmup();
    correctness_spot();

    println!("\n== field multiplication ==");
    bench_field_mul_mine(mul_iters, runs);
    bench_field_mul_k256(mul_iters, runs);

    println!("\n== single inversion ==");
    bench_single_inv_variant("mine_inv_baseline", field_inv, inv_iters, runs);
    bench_single_inv_variant("mine_inv_window4", inv_window4_v1, inv_iters, runs);
    bench_single_inv_k256(inv_iters, runs);

    println!("\n== real normalization workload ==");
    println!("points per run: {}", norm_points);
    bench_norm_variant(NormVariant::SingleBaseline, norm_points, runs);
    bench_norm_variant(NormVariant::SingleWindow4, norm_points, runs);
    bench_norm_variant(NormVariant::BatchAllBaseline, norm_points, runs);
    bench_norm_variant(NormVariant::BatchAllWindow4, norm_points, runs);
    bench_norm_k256(norm_points, runs);

    scaling_sweep();

    println!("\n== timing sanity ==");
    timing_sanity(field_inv, "baseline_inv");
    timing_sanity(inv_window4_v1, "window4_inv");

    println!("\n== phase 2 complete ==");
}
