use gaussian_crypto::{field_inv, field_mul, field_square, P, U256};
use k256::{FieldElement, ProjectivePoint};
use k256::elliptic_curve::BatchNormalize;
use k256::elliptic_curve::sec1::ToEncodedPoint;
use std::hint::black_box;
use std::time::Instant;

const ONE: U256 = [1, 0, 0, 0];
const P_MINUS_2: U256 = [
    0xFFFFFFFEFFFFFC2D,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
];

#[derive(Clone, Copy)]
enum InvChoice {
    Baseline,
    Window4,
    Window5,
}

impl InvChoice {
    fn name(&self) -> &'static str {
        match self {
            InvChoice::Baseline => "baseline",
            InvChoice::Window4 => "window4",
            InvChoice::Window5 => "window5",
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

fn inv_window<const W: usize>(a: U256) -> U256 {
    if a == [0, 0, 0, 0] {
        return [0, 0, 0, 0];
    }

    let table_len = 1usize << W;
    let mut table = vec![[0u64; 4]; table_len];
    table[0] = ONE;
    table[1] = a;
    for i in 2..table_len {
        table[i] = field_mul(table[i - 1], a);
    }

    let exp = u256_to_be32(P_MINUS_2);
    let mut result = ONE;

    let mut bitbuf: u32 = 0;
    let mut bits: usize = 0;
    for byte in exp {
        bitbuf = (bitbuf << 8) | byte as u32;
        bits += 8;
        while bits >= W {
            let shift = bits - W;
            let idx = ((bitbuf >> shift) & ((1u32 << W) - 1)) as usize;
            for _ in 0..W {
                result = field_square(result);
            }
            if idx != 0 {
                result = field_mul(result, table[idx]);
            }
            bitbuf &= (1u32 << shift).wrapping_sub(1);
            bits -= W;
        }
    }

    if bits > 0 {
        for _ in 0..bits {
            result = field_square(result);
        }
        let idx = bitbuf as usize;
        if idx != 0 {
            result = field_mul(result, table[idx]);
        }
    }

    result
}

fn inv_window4(a: U256) -> U256 {
    inv_window::<4>(a)
}

fn inv_window5(a: U256) -> U256 {
    inv_window::<5>(a)
}

fn inv_with(choice: InvChoice, a: U256) -> U256 {
    match choice {
        InvChoice::Baseline => field_inv(a),
        InvChoice::Window4 => inv_window4(a),
        InvChoice::Window5 => inv_window5(a),
    }
}

fn batch_invert_with(values: &[U256], choice: InvChoice) -> Vec<U256> {
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

    let mut inv_acc = inv_with(choice, acc);
    let mut out = vec![[0u64; 4]; n];

    for i in (0..n).rev() {
        out[i] = field_mul(inv_acc, prefix[i]);
        inv_acc = field_mul(inv_acc, values[i]);
    }

    out
}

fn correctness_spot() {
    let mut seed = 0xDEADBEEF12345678u64;
    for _ in 0..3000 {
        let a = rand_nonzero_below_p(&mut seed);
        let inv0 = field_inv(a);
        let inv4 = inv_window4(a);
        let inv5 = inv_window5(a);
        assert_eq!(field_mul(a, inv0), ONE, "baseline inverse failed");
        assert_eq!(field_mul(a, inv4), ONE, "window4 inverse failed");
        assert_eq!(field_mul(a, inv5), ONE, "window5 inverse failed");
    }
    println!("correctness_spot: PASS");
}

fn warmup() {
    let mut seed = 0x5555AAAACCCC9999u64;
    let mut sink = 0u64;
    for _ in 0..20000 {
        let a = rand_nonzero_below_p(&mut seed);
        let b = rand_nonzero_below_p(&mut seed);
        sink ^= field_mul(a, b)[0];
        sink ^= inv_window4(a)[0];
        sink ^= inv_window5(a)[0];
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
            checksum ^= field_mul(black_box(a), black_box(b))[0];
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

fn bench_inv_variant(label: &str, choice: InvChoice, iters: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0x9999AAAABBBBCCCCu64 ^ run as u64;
        let start = Instant::now();
        let mut checksum = 0u64;
        for _ in 0..iters {
            let a = rand_nonzero_below_p(&mut seed);
            checksum ^= inv_with(choice, black_box(a))[0];
        }
        let dur = start.elapsed();
        println!("{} run {:02}: {:?}, checksum={}", label, run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize(&format!("{} summary_sec", label), &results);
}

fn bench_inv_k256(iters: usize, runs: usize) {
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

fn bench_norm_variant(choice: InvChoice, points: usize, runs: usize) {
    let mut results = Vec::with_capacity(runs);
    for run in 0..runs {
        let mut seed = 0xABCD1234EF998877u64 ^ run as u64;
        let xs: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();
        let ys: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();
        let zs: Vec<U256> = (0..points).map(|_| rand_nonzero_below_p(&mut seed)).collect();

        let start = Instant::now();
        let zinv = batch_invert_with(&zs, choice);
        let mut checksum = 0u64;

        for i in 0..points {
            let z2 = field_square(zinv[i]);
            let z3 = field_mul(z2, zinv[i]);
            let ax = field_mul(xs[i], z2);
            let ay = field_mul(ys[i], z3);
            checksum ^= black_box(ax[0] ^ ay[0]);
        }

        let dur = start.elapsed();
        println!("batch_all_{}_norm run {:02}: {:?}, checksum={}", choice.name(), run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }
    summarize(&format!("batch_all_{}_norm summary_sec", choice.name()), &results);
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
        bench_norm_variant(InvChoice::Baseline, points, runs);
        bench_norm_variant(InvChoice::Window4, points, runs);
        bench_norm_variant(InvChoice::Window5, points, runs);
        bench_norm_k256(points, runs);
        println!();
    }
}

fn main() {
    let runs = 20usize;
    let mul_iters = 100_000usize;
    let inv_iters = 10_000usize;
    let norm_points = 8_192usize;

    println!("== final product tournament ==");
    warmup();
    correctness_spot();

    println!("\n== field multiplication ==");
    bench_field_mul_mine(mul_iters, runs);
    bench_field_mul_k256(mul_iters, runs);

    println!("\n== single inversion ==");
    bench_inv_variant("mine_inv_baseline", InvChoice::Baseline, inv_iters, runs);
    bench_inv_variant("mine_inv_window4", InvChoice::Window4, inv_iters, runs);
    bench_inv_variant("mine_inv_window5", InvChoice::Window5, inv_iters, runs);
    bench_inv_k256(inv_iters, runs);

    println!("\n== real normalization workload ==");
    println!("points per run: {}", norm_points);
    bench_norm_variant(InvChoice::Baseline, norm_points, runs);
    bench_norm_variant(InvChoice::Window4, norm_points, runs);
    bench_norm_variant(InvChoice::Window5, norm_points, runs);
    bench_norm_k256(norm_points, runs);

    scaling_sweep();

    println!("\n== final product tournament complete ==");
}
