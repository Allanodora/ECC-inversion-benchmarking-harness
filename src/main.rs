use gaussian_crypto::{field_inv, field_mul, rand_nonzero_below_p, U256};
use std::hint::black_box;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    let mut seed = 0x123456789ABCDEF0u64;
    let mut checksum = 0u64;

    for _ in 0..200 {
        let a: U256 = rand_nonzero_below_p(&mut seed);
        let inv = field_inv(black_box(a));
        let one = field_mul(a, inv);
        checksum ^= black_box(one[0]);
    }

    let duration = start.elapsed();
    println!("Time for 200 correct field inverses + verify mul: {:?}", duration);
    println!("Checksum low limb: {}", checksum);
}
