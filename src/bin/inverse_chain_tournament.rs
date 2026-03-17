use core::hint::black_box;
use std::time::Instant;

use k256::elliptic_curve::ff::{Field, PrimeField};
use k256::FieldElement;

type U256 = [u64; 4];

#[inline(always)]
fn checksum_fe(x: &FieldElement) -> u64 {
    let b: [u8; 32] = x.to_repr().into();
    let w0 = u64::from_be_bytes(b[0..8].try_into().unwrap());
    let w1 = u64::from_be_bytes(b[8..16].try_into().unwrap());
    let w2 = u64::from_be_bytes(b[16..24].try_into().unwrap());
    let w3 = u64::from_be_bytes(b[24..32].try_into().unwrap());
    w0.wrapping_mul(0x9E3779B97F4A7C15) ^ w1.rotate_left(13) ^ w2.rotate_left(29) ^ w3.rotate_left(47)
}

#[inline(always)]
fn fe_equal(a: &FieldElement, b: &FieldElement) -> bool {
    let ab: [u8; 32] = a.to_repr().into();
    let bb: [u8; 32] = b.to_repr().into();
    ab == bb
}

#[inline(always)]
fn xorshift64star(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    x.wrapping_mul(0x2545F4914F6CDD1D)
}

#[inline(always)]
fn sqn(mut x: FieldElement, n: usize) -> FieldElement {
    for _ in 0..n {
        x = x.square();
    }
    x
}

fn make_inputs(n: usize) -> Vec<FieldElement> {
    let mut s = 0x1234_5678_9ABC_DEF0u64;
    let mut out = Vec::with_capacity(n);

    while out.len() < n {
        let mut bytes = [0u8; 32];
        for chunk in bytes.chunks_exact_mut(8) {
            chunk.copy_from_slice(&xorshift64star(&mut s).to_be_bytes());
        }

        if let Some(fe) = Option::<FieldElement>::from(FieldElement::from_repr(bytes.into())) {
            if !bool::from(fe.is_zero()) {
                out.push(fe);
            }
        }
    }

    out
}

#[inline(always)]
fn k256_inv_fe(x: &FieldElement) -> FieldElement {
    x.invert().unwrap()
}

/*
Translated fixed addition chain for secp256k1 field inversion.

This follows the classic secp256k1_fe_inv chain shape:
[1], [2], 3, 6, 9, 11, [22], 44, 88, 176, 220, [223]
then final assembly with:
+23 squares, *x22, +5 squares, *a, +3 squares, *x2, +2 squares, *a
*/
#[inline(always)]
fn fixed_chain_inv_fe(a: &FieldElement) -> FieldElement {
    let x2 = a.square() * *a;               // a^(2^2 - 1)  = a^3
    let x3 = x2.square() * *a;              // a^(2^3 - 1)  = a^7

    let x6 = sqn(x3, 3) * x3;               // a^(2^6 - 1)
    let x9 = sqn(x6, 3) * x3;               // a^(2^9 - 1)
    let x11 = sqn(x9, 2) * x2;              // a^(2^11 - 1)
    let x22 = sqn(x11, 11) * x11;           // a^(2^22 - 1)
    let x44 = sqn(x22, 22) * x22;           // a^(2^44 - 1)
    let x88 = sqn(x44, 44) * x44;           // a^(2^88 - 1)
    let x176 = sqn(x88, 88) * x88;          // a^(2^176 - 1)
    let x220 = sqn(x176, 44) * x44;         // a^(2^220 - 1)
    let x223 = sqn(x220, 3) * x3;           // a^(2^223 - 1)

    let mut t1 = sqn(x223, 23) * x22;
    t1 = sqn(t1, 5) * *a;
    t1 = sqn(t1, 3) * x2;
    sqn(t1, 2) * *a
}

#[derive(Clone, Copy, Debug)]
struct Stats {
    min: f64,
    median: f64,
    mean: f64,
    max: f64,
}

fn summarize(times: &mut [f64]) -> Stats {
    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min = times[0];
    let max = times[times.len() - 1];
    let median = if times.len() % 2 == 0 {
        let i = times.len() / 2;
        (times[i - 1] + times[i]) * 0.5
    } else {
        times[times.len() / 2]
    };
    let mean = times.iter().sum::<f64>() / times.len() as f64;
    Stats { min, median, mean, max }
}

fn bench_single<F>(name: &str, runs: usize, inputs: &[FieldElement], mut f: F)
where
    F: FnMut(&FieldElement) -> FieldElement,
{
    let mut warm = 0u64;
    for x in inputs.iter().take(64) {
        let y = black_box(f(black_box(x)));
        warm ^= checksum_fe(&y);
    }
    black_box(warm);

    let mut times = Vec::with_capacity(runs);

    for run in 0..runs {
        let t0 = Instant::now();
        let mut chk = 0u64;

        for x in inputs {
            let y = black_box(f(black_box(x)));
            chk ^= checksum_fe(&y);
        }

        let dt = t0.elapsed().as_secs_f64();
        println!(
            "{} run {:02}: {:>10.6}ms checksum={}",
            name,
            run + 1,
            dt * 1000.0,
            chk
        );
        times.push(dt);
    }

    let s = summarize(&mut times);
    println!(
        "{} summary_sec min={:.9} median={:.9} mean={:.9} max={:.9}",
        name, s.min, s.median, s.mean, s.max
    );
}

fn correctness_spot(inputs: &[FieldElement]) {
    for (i, x) in inputs.iter().take(64).enumerate() {
        let a = fixed_chain_inv_fe(x);
        let b = k256_inv_fe(x);
        assert!(fe_equal(&a, &b), "fixed-chain mismatch at {}", i);
    }
    println!("correctness_spot: PASS");
}

fn main() {
    const RUNS: usize = 20;
    const N_SINGLE: usize = 4096;

    println!("== inverse chain tournament ==");

    let inputs = make_inputs(N_SINGLE);

    correctness_spot(&inputs);

    println!("\n== standalone inverse ==");

    bench_single(
        "fixed_chain_inv",
        RUNS,
        &inputs,
        |x| fixed_chain_inv_fe(x),
    );

    bench_single(
        "k256_inv",
        RUNS,
        &inputs,
        |x| k256_inv_fe(x),
    );

    println!("\n== inverse chain tournament complete ==");
}
