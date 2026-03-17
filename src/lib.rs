use num_bigint::BigUint;
use num_traits::{One, Zero};

pub type U256 = [u64; 4];
pub type U512 = [u64; 8];

pub const P: U256 = [
    0xFFFFFFFEFFFFFC2F,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
];

const MASK64: u128 = 0xFFFF_FFFF_FFFF_FFFF;
const C977: u128 = 977;
const CSH: u128 = 1u128 << 32;

#[inline]
pub fn cmp_u256(a: &U256, b: &U256) -> core::cmp::Ordering {
    for i in (0..4).rev() {
        if a[i] < b[i] {
            return core::cmp::Ordering::Less;
        } else if a[i] > b[i] {
            return core::cmp::Ordering::Greater;
        }
    }
    core::cmp::Ordering::Equal
}

#[inline]
pub fn add_u256(a: U256, b: U256) -> (U256, u64) {
    let mut out = [0u64; 4];
    let mut carry = 0u64;
    for i in 0..4 {
        let (r1, c1) = a[i].overflowing_add(b[i]);
        let (r2, c2) = r1.overflowing_add(carry);
        out[i] = r2;
        carry = (c1 as u64) | (c2 as u64);
    }
    (out, carry)
}

#[inline]
pub fn sub_u256(a: U256, b: U256) -> (U256, u64) {
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

#[inline]
pub fn mul_u256_wide(a: U256, b: U256) -> U512 {
    let mut out = [0u64; 8];

    for i in 0..4 {
        let mut carry = 0u128;
        for j in 0..4 {
            let k = i + j;
            let t = (a[i] as u128)
                .wrapping_mul(b[j] as u128)
                .wrapping_add(out[k] as u128)
                .wrapping_add(carry);
            out[k] = t as u64;
            carry = t >> 64;
        }

        let mut k = i + 4;
        let mut c = carry;
        while c != 0 && k < 8 {
            let t = (out[k] as u128).wrapping_add(c);
            out[k] = t as u64;
            c = t >> 64;
            k += 1;
        }
    }

    out
}

#[inline]
pub fn lo_u256(x: U512) -> U256 {
    [x[0], x[1], x[2], x[3]]
}

#[inline]
pub fn u256_to_big(x: U256) -> BigUint {
    let mut bytes = Vec::with_capacity(32);
    for limb in x {
        bytes.extend_from_slice(&limb.to_le_bytes());
    }
    BigUint::from_bytes_le(&bytes)
}

#[inline]
pub fn u512_to_big(x: U512) -> BigUint {
    let mut bytes = Vec::with_capacity(64);
    for limb in x {
        bytes.extend_from_slice(&limb.to_le_bytes());
    }
    BigUint::from_bytes_le(&bytes)
}

#[inline]
pub fn big_to_u256(x: &BigUint) -> U256 {
    let bytes = x.to_bytes_le();
    let mut out = [0u64; 4];
    for i in 0..4 {
        let mut limb = [0u8; 8];
        for j in 0..8 {
            let idx = i * 8 + j;
            if idx < bytes.len() {
                limb[j] = bytes[idx];
            }
        }
        out[i] = u64::from_le_bytes(limb);
    }
    out
}

#[inline]
pub fn p_big() -> BigUint {
    u256_to_big(P)
}

#[inline]
pub fn reduce_once(x: U256) -> U256 {
    if cmp_u256(&x, &P) != core::cmp::Ordering::Less {
        let (y, _) = sub_u256(x, P);
        y
    } else {
        x
    }
}

#[inline]
fn normalize_8(acc: &mut [u128; 8]) {
    for i in 0..7 {
        let carry = acc[i] >> 64;
        acc[i] &= MASK64;
        acc[i + 1] = acc[i + 1].wrapping_add(carry);
    }
    acc[7] &= MASK64;
}

#[inline]
pub fn mod_p_u512(x: U512) -> U256 {
    let mut acc = [0u128; 8];
    for i in 0..8 {
        acc[i] = x[i] as u128;
    }

    normalize_8(&mut acc);

    loop {
        let mut changed = false;

        for k in 4..8 {
            let v = acc[k];
            if v != 0 {
                acc[k] = 0;
                let base = k - 4;
                acc[base] = acc[base].wrapping_add(v.wrapping_mul(C977));
                acc[base] = acc[base].wrapping_add(v.wrapping_mul(CSH));
                changed = true;
            }
        }

        normalize_8(&mut acc);

        if !changed {
            break;
        }
    }

    let mut out = [acc[0] as u64, acc[1] as u64, acc[2] as u64, acc[3] as u64];
    while cmp_u256(&out, &P) != core::cmp::Ordering::Less {
        let (t, _) = sub_u256(out, P);
        out = t;
    }
    out
}

#[inline]
pub fn field_mul(a: U256, b: U256) -> U256 {
    mod_p_u512(mul_u256_wide(a, b))
}

#[inline]
pub fn field_square(a: U256) -> U256 {
    field_mul(a, a)
}

#[inline]
pub fn field_pow(base: U256, exp: &BigUint) -> U256 {
    let mut result = [1u64, 0, 0, 0];
    let mut cur = base;
    let mut e = exp.clone();

    while !e.is_zero() {
        if (&e & BigUint::one()) == BigUint::one() {
            result = field_mul(result, cur);
        }
        cur = field_square(cur);
        e >>= 1usize;
    }

    result
}

#[inline]
pub fn field_inv(a: U256) -> U256 {
    if a == [0, 0, 0, 0] {
        return [0, 0, 0, 0];
    }
    let exp = p_big() - BigUint::from(2u32);
    field_pow(a, &exp)
}

#[inline]
pub fn lcg(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *seed
}

#[inline]
pub fn rand_u256(seed: &mut u64) -> U256 {
    [lcg(seed), lcg(seed), lcg(seed), lcg(seed)]
}

#[inline]
pub fn rand_nonzero_below_p(seed: &mut u64) -> U256 {
    loop {
        let x = reduce_once(rand_u256(seed));
        if x != [0, 0, 0, 0] {
            return x;
        }
    }
}
