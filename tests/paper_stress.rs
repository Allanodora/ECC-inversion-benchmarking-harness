use gaussian_crypto::{field_inv, field_mul, P, U256};

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

#[test]
fn stress_field_mul_1m() {
    let mut seed = 0x1111222233334444u64;
    let mut checksum = 0u64;

    for _ in 0..1_000_000u64 {
        let a = rand_nonzero_below_p(&mut seed);
        let b = rand_nonzero_below_p(&mut seed);
        let c = field_mul(a, b);
        checksum ^= c[0];
    }

    println!("stress_field_mul_1m checksum={}", checksum);
}

#[test]
fn stress_field_inv_100k() {
    let mut seed = 0x5555666677778888u64;
    let mut checksum = 0u64;

    for _ in 0..100_000u64 {
        let a = rand_nonzero_below_p(&mut seed);
        let inv = field_inv(a);
        let one = field_mul(a, inv);
        assert_eq!(one, [1,0,0,0]);
        checksum ^= inv[0];
    }

    println!("stress_field_inv_100k checksum={}", checksum);
}

#[test]
fn worst_case_vectors() {
    let cases: [U256; 10] = [
        [0xFFFFFFFFFFFFFFFF; 4],
        [0xAAAAAAAAAAAAAAAA; 4],
        [0x5555555555555555; 4],
        [P[0].wrapping_sub(1), P[1], P[2], P[3]],
        [P[0].wrapping_sub(2), P[1], P[2], P[3]],
        [P[0].wrapping_sub(3), P[1], P[2], P[3]],
        [u64::MAX, u64::MAX, 0, 0],
        [u64::MAX, 0, u64::MAX, 0],
        [0, u64::MAX, 0, u64::MAX],
        [1, 0, 0, 0],
    ];

    for &a in &cases {
        for &b in &cases {
            let aa = reduce_once(a);
            let bb = reduce_once(b);
            let c = field_mul(aa, bb);

            if aa != [0,0,0,0] {
                let inv = field_inv(aa);
                let one = field_mul(aa, inv);
                assert_eq!(one, [1,0,0,0]);
            }

            let _ = c;
        }
    }
}
