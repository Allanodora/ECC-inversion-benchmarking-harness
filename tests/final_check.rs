use gaussian_crypto::{
    add_u256, big_to_u256, cmp_u256, field_inv, field_mul, lo_u256, mod_p_u512,
    mul_u256_wide, p_big, rand_nonzero_below_p, rand_u256, reduce_once, sub_u256,
    u256_to_big, u512_to_big, P, U256,
};
use num_bigint::BigUint;
use num_traits::{One, Zero};

fn two_256_mask() -> BigUint {
    (BigUint::one() << 256usize) - BigUint::one()
}

#[test]
fn cmp_edges() {
    assert!(cmp_u256(&[0,0,0,0], &[0,0,0,0]).is_eq());
    assert!(cmp_u256(&[1,0,0,0], &[2,0,0,0]).is_lt());
    assert!(cmp_u256(&[0,0,0,2], &[0,0,0,1]).is_gt());
}

#[test]
fn add_sub_random_100k() {
    let mut seed = 1u64;
    let mask = two_256_mask();

    for _ in 0..100_000 {
        let a = rand_u256(&mut seed);
        let b = rand_u256(&mut seed);

        let (sum, carry) = add_u256(a, b);
        let sum_big = u256_to_big(a) + u256_to_big(b);
        let expected_lo = &sum_big & &mask;
        let expected_carry = if (sum_big >> 256usize).is_zero() { 0 } else { 1 };
        assert_eq!(sum, big_to_u256(&expected_lo));
        assert_eq!(carry, expected_carry);

        let (diff, borrow) = sub_u256(a, b);
        let a_big = u256_to_big(a);
        let b_big = u256_to_big(b);
        if a_big >= b_big {
            assert_eq!(borrow, 0);
            assert_eq!(diff, big_to_u256(&(a_big - b_big)));
        } else {
            assert_eq!(borrow, 1);
        }
    }
}

#[test]
fn wide_mul_random_100k() {
    let mut seed = 7u64;

    for _ in 0..100_000 {
        let a = rand_u256(&mut seed);
        let b = rand_u256(&mut seed);
        let got = u512_to_big(mul_u256_wide(a, b));
        let expected = u256_to_big(a) * u256_to_big(b);
        assert_eq!(got, expected);
    }
}

#[test]
fn low_half_random_100k() {
    let mut seed = 123u64;
    let mask = two_256_mask();

    for _ in 0..100_000 {
        let a = rand_u256(&mut seed);
        let b = rand_u256(&mut seed);
        let wide = mul_u256_wide(a, b);
        let low = lo_u256(wide);
        let expected_low = (u256_to_big(a) * u256_to_big(b)) & &mask;
        assert_eq!(low, big_to_u256(&expected_low));
    }
}

#[test]
fn mod_p_random_50k() {
    let mut seed = 0x7777777777777777u64;
    let p = p_big();

    for _ in 0..50_000 {
        let a = rand_u256(&mut seed);
        let b = rand_u256(&mut seed);
        let wide = mul_u256_wide(a, b);
        let got = mod_p_u512(wide);
        let expected = (u256_to_big(a) * u256_to_big(b)) % &p;
        assert_eq!(got, big_to_u256(&expected));
    }
}

#[test]
fn field_mul_random_50k() {
    let mut seed = 0xABCDEFFF01234567u64;
    let p = p_big();

    for _ in 0..50_000 {
        let a = reduce_once(rand_u256(&mut seed));
        let b = reduce_once(rand_u256(&mut seed));
        let got = field_mul(a, b);
        let expected = (u256_to_big(a) * u256_to_big(b)) % &p;
        assert_eq!(got, big_to_u256(&expected));
    }
}

#[test]
fn field_inv_random_1000() {
    let mut seed = 0xDEADBEEFCAFEBABEu64;

    for _ in 0..1000 {
        let a: U256 = rand_nonzero_below_p(&mut seed);
        let inv = field_inv(a);
        let one = field_mul(a, inv);
        assert_eq!(one, [1,0,0,0], "inverse failed for {:?}", a);
    }
}

#[test]
fn reduce_edges() {
    assert_eq!(reduce_once([0,0,0,0]), [0,0,0,0]);
    assert_eq!(reduce_once(P), [0,0,0,0]);

    let (p_plus_1, carry) = add_u256(P, [1,0,0,0]);
    assert_eq!(carry, 0);
    assert_eq!(reduce_once(p_plus_1), [1,0,0,0]);
}
