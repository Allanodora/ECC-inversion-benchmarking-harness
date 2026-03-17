use core::hint::black_box;
use std::time::Instant;

use k256::FieldElement;
use k256::elliptic_curve::ff::PrimeField;

type U256 = [u64;4];

#[inline(always)]
fn checksum_u256(x:&U256)->u64{
    x[0].wrapping_mul(0x9E3779B97F4A7C15)
        ^ x[1].rotate_left(13)
        ^ x[2].rotate_left(29)
        ^ x[3].rotate_left(47)
}

#[inline(always)]
fn cmp_u256(a:&U256,b:&U256)->core::cmp::Ordering{
    if a[3]!=b[3]{return a[3].cmp(&b[3]);}
    if a[2]!=b[2]{return a[2].cmp(&b[2]);}
    if a[1]!=b[1]{return a[1].cmp(&b[1]);}
    a[0].cmp(&b[0])
}

/* ===============================
   U256 ↔ FieldElement conversion
   =============================== */

#[inline(always)]
fn u256_to_bytes(x:U256)->[u8;32]{
    let mut out=[0u8;32];
    out[0..8].copy_from_slice(&x[3].to_be_bytes());
    out[8..16].copy_from_slice(&x[2].to_be_bytes());
    out[16..24].copy_from_slice(&x[1].to_be_bytes());
    out[24..32].copy_from_slice(&x[0].to_be_bytes());
    out
}

#[inline(always)]
fn bytes_to_u256(b:[u8;32])->U256{
    [
        u64::from_be_bytes(b[24..32].try_into().unwrap()),
        u64::from_be_bytes(b[16..24].try_into().unwrap()),
        u64::from_be_bytes(b[8..16].try_into().unwrap()),
        u64::from_be_bytes(b[0..8].try_into().unwrap())
    ]
}

#[inline(always)]
fn fe_from_u256(x:U256)->FieldElement{
    let bytes=u256_to_bytes(x);
    FieldElement::from_repr(bytes.into()).unwrap()
}

#[inline(always)]
fn u256_from_fe(x:FieldElement)->U256{
    let bytes: [u8;32] = x.to_repr().into();
    bytes_to_u256(bytes)
}

#[inline(always)]
fn k256_inv_u256(x:&U256)->U256{
    let fe=fe_from_u256(*x);
    let inv=fe.invert().unwrap();
    u256_from_fe(inv)
}

/* ===============================
   Placeholder algorithms
   =============================== */

#[inline(always)]
fn mine_inv_window5(x:&U256)->U256{
    k256_inv_u256(x)
}

#[inline(always)]
fn mine_inv_fixed_chain(x:&U256)->U256{
    mine_inv_window5(x)
}

fn mine_batch_invert_all(xs:&[U256])->Vec<U256>{
    xs.iter().map(|x|mine_inv_window5(x)).collect()
}

/* ===============================
   Random generator
   =============================== */

#[inline(always)]
fn xorshift64star(state:&mut u64)->u64{
    let mut x=*state;
    x^=x>>12;
    x^=x<<25;
    x^=x>>27;
    *state=x;
    x.wrapping_mul(0x2545F4914F6CDD1D)
}

fn make_inputs(n:usize)->Vec<U256>{
    let mut s=0x123456789ABCDEF0u64;
    let mut out=Vec::with_capacity(n);

    for _ in 0..n{
        let mut x=[
            xorshift64star(&mut s),
            xorshift64star(&mut s),
            xorshift64star(&mut s),
            xorshift64star(&mut s)
        ];

        if x==[0,0,0,0]{x[0]=1;}
        x[3]&=0x7FFFFFFFFFFFFFFF;

        out.push(x);
    }

    out
}

/* ===============================
   Benchmark helpers
   =============================== */

#[derive(Clone,Copy,Debug)]
struct Stats{
    min:f64,
    median:f64,
    mean:f64,
    max:f64
}

fn summarize(times:&mut[f64])->Stats{
    times.sort_by(|a,b|a.partial_cmp(b).unwrap());

    let min=times[0];
    let max=times[times.len()-1];

    let median=if times.len()%2==0{
        let i=times.len()/2;
        (times[i-1]+times[i])*0.5
    }else{
        times[times.len()/2]
    };

    let mean=times.iter().sum::<f64>()/times.len() as f64;

    Stats{min,median,mean,max}
}

fn bench_single<F>(name:&str,runs:usize,inputs:&[U256],mut f:F)
where F:FnMut(&U256)->U256{

    let mut warm=0u64;
    for x in inputs.iter().take(64){
        let y=black_box(f(black_box(x)));
        warm^=checksum_u256(&y);
    }
    black_box(warm);

    let mut times=Vec::with_capacity(runs);

    for run in 0..runs{
        let t0=Instant::now();
        let mut chk=0u64;

        for x in inputs{
            let y=black_box(f(black_box(x)));
            chk^=checksum_u256(&y);
        }

        let dt=t0.elapsed().as_secs_f64();

        println!(
            "{} run {:02}: {:>10.6}ms checksum={}",
            name,
            run+1,
            dt*1000.0,
            chk
        );

        times.push(dt);
    }

    let s=summarize(&mut times);

    println!(
        "{} summary_sec min={:.9} median={:.9} mean={:.9} max={:.9}",
        name,s.min,s.median,s.mean,s.max
    );
}

fn main(){

    const RUNS:usize=20;
    const N_SINGLE:usize=4096;

    println!("== inverse only tournament v2 ==");

    let inputs=make_inputs(N_SINGLE);

    bench_single(
        "mine_inv_window5",
        RUNS,
        &inputs,
        |x|mine_inv_window5(x)
    );

    bench_single(
        "mine_inv_fixed_chain",
        RUNS,
        &inputs,
        |x|mine_inv_fixed_chain(x)
    );

    bench_single(
        "k256_inv",
        RUNS,
        &inputs,
        |x|k256_inv_u256(x)
    );

    println!("\n== inverse only tournament v2 complete ==");
}
