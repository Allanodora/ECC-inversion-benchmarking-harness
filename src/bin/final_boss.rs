use k256::ecdsa::{
    signature::{Signer as KSigner, Verifier as KVerifier},
    Signature as KSignature, SigningKey as KSigningKey, VerifyingKey as KVerifyingKey,
};
use secp256k1::{ecdsa::Signature as SSignature, Message, PublicKey, Secp256k1, SecretKey};
use sha2::{Digest, Sha256};
use std::hint::black_box;
use std::time::Instant;

fn hash32(msg: &[u8]) -> [u8; 32] {
    let mut hasher = Sha256::new();
    hasher.update(msg);
    let out = hasher.finalize();
    let mut arr = [0u8; 32];
    arr.copy_from_slice(&out);
    arr
}

fn deterministic_secret(i: u64) -> [u8; 32] {
    let mut seed = [0u8; 32];
    seed[0..8].copy_from_slice(&i.to_be_bytes());
    let mut out = hash32(&seed);
    if out.iter().all(|&b| b == 0) {
        out[31] = 1;
    }
    out
}

fn make_msg(run: usize, i: usize, sk_bytes: [u8; 32]) -> [u8; 32] {
    let mut data = Vec::with_capacity(40);
    data.extend_from_slice(&(((run as u64) << 32) | i as u64).to_be_bytes());
    data.extend_from_slice(&sk_bytes);
    hash32(&data)
}

fn bench_k256_keygen(iters: usize) {
    let runs = 6usize;
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let start = Instant::now();
        let mut checksum = 0u64;

        for i in 0..iters {
            let sk_bytes = deterministic_secret(((run as u64) << 32) | i as u64);
            let signing_key = KSigningKey::from_bytes((&sk_bytes).into()).expect("valid k256 sk");
            let verify_key: KVerifyingKey = *signing_key.verifying_key();
            let enc = verify_key.to_encoded_point(true);
            let bytes = enc.as_bytes();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[bytes.len() - 8..]);
            checksum ^= u64::from_be_bytes(low);
            black_box(verify_key);
        }

        let dur = start.elapsed();
        println!("k256 keygen run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("k256 keygen summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_secp_keygen(iters: usize) {
    let runs = 6usize;
    let secp = Secp256k1::new();
    let mut results = Vec::with_capacity(runs);

    for run in 0..runs {
        let start = Instant::now();
        let mut checksum = 0u64;

        for i in 0..iters {
            let sk_bytes = deterministic_secret(((run as u64) << 32) | i as u64);
            let sk = SecretKey::from_byte_array(sk_bytes).expect("valid secp sk");
            let pk = PublicKey::from_secret_key_global(&sk);
            let ser = pk.serialize();
            let mut low = [0u8; 8];
            low.copy_from_slice(&ser[ser.len() - 8..]);
            checksum ^= u64::from_be_bytes(low);
            black_box((&secp, pk));
        }

        let dur = start.elapsed();
        println!("secp256k1 keygen run {}: {:?}, checksum={}", run + 1, dur, checksum);
        results.push(dur.as_secs_f64());
    }

    let min = results.iter().copied().fold(f64::INFINITY, f64::min);
    let max = results.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = results.iter().sum::<f64>() / results.len() as f64;
    println!("secp256k1 keygen summary_sec min={:.9} mean={:.9} max={:.9}", min, mean, max);
}

fn bench_k256_sign_verify(iters: usize) {
    let runs = 5usize;
    let mut results_sign = Vec::with_capacity(runs);
    let mut results_verify = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut sigs = Vec::with_capacity(iters);
        let mut vks = Vec::with_capacity(iters);
        let mut msgs = Vec::with_capacity(iters);

        for i in 0..iters {
            let sk_bytes = deterministic_secret(((run as u64) << 32) | i as u64);
            let signing_key = KSigningKey::from_bytes((&sk_bytes).into()).expect("valid k256 sk");
            let verify_key: KVerifyingKey = *signing_key.verifying_key();
            let msg = make_msg(run, i, sk_bytes);

            let start = Instant::now();
            let sig: KSignature = signing_key.sign(&msg);
            let dur = start.elapsed();

            sigs.push(sig);
            vks.push(verify_key);
            msgs.push(msg);
            results_sign.push(dur.as_secs_f64());
        }

        let start_v = Instant::now();
        let mut checksum = 0u64;
        for i in 0..iters {
            vks[i].verify(&msgs[i], &sigs[i]).expect("k256 verify");
            let bytes = sigs[i].to_bytes();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[bytes.len() - 8..]);
            checksum ^= u64::from_be_bytes(low);
        }
        let dur_v = start_v.elapsed();
        println!("k256 sign/verify run {}: verify {:?}, checksum={}", run + 1, dur_v, checksum);
        results_verify.push(dur_v.as_secs_f64());
    }

    let sign_mean = results_sign.iter().sum::<f64>() / results_sign.len() as f64;
    let verify_mean = results_verify.iter().sum::<f64>() / results_verify.len() as f64;
    println!("k256 sign mean_sec_per_op={:.9}", sign_mean);
    println!("k256 verify mean_sec_per_run={:.9}", verify_mean);
}

fn bench_secp_sign_verify(iters: usize) {
    let runs = 5usize;
    let secp = Secp256k1::new();
    let mut results_sign = Vec::with_capacity(runs);
    let mut results_verify = Vec::with_capacity(runs);

    for run in 0..runs {
        let mut sigs = Vec::with_capacity(iters);
        let mut pks = Vec::with_capacity(iters);
        let mut msgs = Vec::with_capacity(iters);

        for i in 0..iters {
            let sk_bytes = deterministic_secret(((run as u64) << 32) | i as u64);
            let sk = SecretKey::from_byte_array(sk_bytes).expect("valid secp sk");
            let pk = PublicKey::from_secret_key_global(&sk);
            let digest = make_msg(run, i, sk_bytes);
            let msg = Message::from_digest(digest);

            let start = Instant::now();
            let sig: SSignature = secp.sign_ecdsa(msg, &sk);
            let dur = start.elapsed();

            sigs.push(sig);
            pks.push(pk);
            msgs.push(msg);
            results_sign.push(dur.as_secs_f64());
        }

        let start_v = Instant::now();
        let mut checksum = 0u64;
        for i in 0..iters {
            secp.verify_ecdsa(msgs[i], &sigs[i], &pks[i]).expect("secp verify");
            let bytes = sigs[i].serialize_compact();
            let mut low = [0u8; 8];
            low.copy_from_slice(&bytes[bytes.len() - 8..]);
            checksum ^= u64::from_be_bytes(low);
        }
        let dur_v = start_v.elapsed();
        println!("secp256k1 sign/verify run {}: verify {:?}, checksum={}", run + 1, dur_v, checksum);
        results_verify.push(dur_v.as_secs_f64());
    }

    let sign_mean = results_sign.iter().sum::<f64>() / results_sign.len() as f64;
    let verify_mean = results_verify.iter().sum::<f64>() / results_verify.len() as f64;
    println!("secp256k1 sign mean_sec_per_op={:.9}", sign_mean);
    println!("secp256k1 verify mean_sec_per_run={:.9}", verify_mean);
}

fn main() {
    println!("== final boss reference tournament ==");
    println!("These are reference baselines for tasks your code does not implement yet.");

    let keygen_iters = 10_000usize;
    let sign_verify_iters = 2_000usize;

    println!("\n== key generation ==");
    bench_k256_keygen(keygen_iters);
    bench_secp_keygen(keygen_iters);

    println!("\n== sign + verify ==");
    bench_k256_sign_verify(sign_verify_iters);
    bench_secp_sign_verify(sign_verify_iters);
}
