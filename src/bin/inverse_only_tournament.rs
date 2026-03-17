use gaussian_crypto::{field_inv, field_mul, field_square, P, U256};

use k256::{FieldElement, ProjectivePoint, Scalar};
use k256::elliptic_curve::BatchNormalize;
use k256::elliptic_curve::sec1::ToEncodedPoint;

use std::hint::black_box;
use std::time::Instant;

const ONE: U256 = [1,0,0,0];

const P_MINUS_2: U256 = [
    0xFFFFFFFEFFFFFC2D,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF,
    0xFFFFFFFFFFFFFFFF
];

#[derive(Clone,Copy)]
enum InvChoice {
    Baseline,
    Window3,
    Window4,
    Window5,
    Window6
}

impl InvChoice {
    fn name(&self)->&'static str{
        match self{
            InvChoice::Baseline=>"baseline",
            InvChoice::Window3=>"window3",
            InvChoice::Window4=>"window4",
            InvChoice::Window5=>"window5",
            InvChoice::Window6=>"window6"
        }
    }
}

fn lcg(seed:&mut u64)->u64{
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *seed
}

fn rand_u256(seed:&mut u64)->U256{
    [lcg(seed),lcg(seed),lcg(seed),lcg(seed)]
}

fn cmp_u256(a:&U256,b:&U256)->core::cmp::Ordering{
    for i in (0..4).rev(){
        if a[i]<b[i]{return core::cmp::Ordering::Less}
        if a[i]>b[i]{return core::cmp::Ordering::Greater}
    }
    core::cmp::Ordering::Equal
}

fn sub_u256(a:U256,b:U256)->(U256,u64){
    let mut out=[0u64;4];
    let mut borrow=0u64;

    for i in 0..4{
        let (r1,b1)=a[i].overflowing_sub(b[i]);
        let (r2,b2)=r1.overflowing_sub(borrow);
        out[i]=r2;
        borrow=(b1 as u64)|(b2 as u64);
    }

    (out,borrow)
}

fn reduce_once(x:U256)->U256{
    if cmp_u256(&x,&P)!=core::cmp::Ordering::Less{
        let (y,_) = sub_u256(x,P);
        y
    }else{
        x
    }
}

fn rand_nonzero_below_p(seed:&mut u64)->U256{
    loop{
        let x=reduce_once(rand_u256(seed));
        if x!=[0,0,0,0]{
            return x
        }
    }
}

fn u256_to_be32(x:U256)->[u8;32]{
    let mut out=[0u8;32];

    for i in 0..4{
        out[(3-i)*8..(4-i)*8].copy_from_slice(&x[i].to_be_bytes());
    }

    out
}

fn to_fe(x:U256)->FieldElement{
    let bytes=u256_to_be32(x);
    let ct=FieldElement::from_bytes((&bytes).into());
    Option::<FieldElement>::from(ct).expect("valid field element")
}

fn median(values:&[f64])->f64{
    let mut v=values.to_vec();
    v.sort_by(|a,b|a.partial_cmp(b).unwrap());

    if v.len()%2==0{
        (v[v.len()/2-1]+v[v.len()/2])/2.0
    }else{
        v[v.len()/2]
    }
}

fn summarize(label:&str,values:&[f64]){
    let mut v=values.to_vec();
    v.sort_by(|a,b|a.partial_cmp(b).unwrap());

    let min=v[0];
    let max=v[v.len()-1];
    let mean=v.iter().sum::<f64>()/v.len() as f64;
    let med=median(values);

    println!(
        "{} min={:.9} median={:.9} mean={:.9} max={:.9}",
        label,min,med,mean,max
    );
}

fn inv_window<const W:usize>(a:U256)->U256{

    if a==[0,0,0,0]{
        return [0,0,0,0]
    }

    let table_len=1usize<<W;

    let mut table=vec![[0u64;4];table_len];

    table[0]=ONE;
    table[1]=a;

    for i in 2..table_len{
        table[i]=field_mul(table[i-1],a);
    }

    let exp=u256_to_be32(P_MINUS_2);

    let mut result=ONE;

    let mut bitbuf: u32=0;
    let mut bits: usize=0;

    for byte in exp{

        bitbuf=(bitbuf<<8)|(byte as u32);
        bits+=8;

        while bits>=W{

            let shift=bits-W;

            let idx=((bitbuf>>shift)&((1u32<<W)-1)) as usize;

            for _ in 0..W{
                result=field_square(result);
            }

            if idx!=0{
                result=field_mul(result,table[idx]);
            }

            bitbuf&=(1u32<<shift).wrapping_sub(1);

            bits-=W;
        }
    }

    if bits>0{

        for _ in 0..bits{
            result=field_square(result);
        }

        let idx=bitbuf as usize;

        if idx!=0{
            result=field_mul(result,table[idx]);
        }
    }

    result
}

fn inv_with(choice:InvChoice,a:U256)->U256{
    match choice{
        InvChoice::Baseline=>field_inv(a),
        InvChoice::Window3=>inv_window::<3>(a),
        InvChoice::Window4=>inv_window::<4>(a),
        InvChoice::Window5=>inv_window::<5>(a),
        InvChoice::Window6=>inv_window::<6>(a)
    }
}

fn correctness_spot(){

    let mut seed=0xDEADBEEF12345678u64;

    let choices=[
        InvChoice::Baseline,
        InvChoice::Window3,
        InvChoice::Window4,
        InvChoice::Window5,
        InvChoice::Window6
    ];

    for _ in 0..3000{

        let a=rand_nonzero_below_p(&mut seed);

        for choice in choices{

            let inv=inv_with(choice,a);

            let one=field_mul(a,inv);

            assert_eq!(one,ONE,"{} inverse failed",choice.name());
        }
    }

    println!("correctness_spot: PASS");
}

fn warmup(){

    let mut seed=0x5555AAAACCCC9999u64;

    let mut sink=0u64;

    for _ in 0..20000{

        let a=rand_nonzero_below_p(&mut seed);

        sink^=inv_with(InvChoice::Window3,a)[0];
        sink^=inv_with(InvChoice::Window4,a)[0];
        sink^=inv_with(InvChoice::Window5,a)[0];
        sink^=inv_with(InvChoice::Window6,a)[0];
    }

    black_box(sink);
}

fn bench_inv_variant(choice:InvChoice,iters:usize,runs:usize){

    let mut results=Vec::with_capacity(runs);

    for run in 0..runs{

        let mut seed=0x9999AAAABBBBCCCCu64 ^ run as u64;

        let start=Instant::now();

        let mut checksum=0u64;

        for _ in 0..iters{

            let a=rand_nonzero_below_p(&mut seed);

            checksum^=inv_with(choice,black_box(a))[0];
        }

        let dur=start.elapsed();

        println!(
            "mine_inv_{} run {:02}: {:?}, checksum={}",
            choice.name(),run+1,dur,checksum
        );

        results.push(dur.as_secs_f64());
    }

    summarize(&format!("mine_inv_{} summary_sec",choice.name()),&results);
}

fn bench_inv_k256(iters:usize,runs:usize){

    let mut results=Vec::with_capacity(runs);

    for run in 0..runs{

        let mut seed=0x9999AAAABBBBCCCCu64 ^ run as u64;

        let start=Instant::now();

        let mut checksum=0u64;

        for _ in 0..iters{

            let a=to_fe(rand_nonzero_below_p(&mut seed));

            let inv_ct=black_box(a).invert();

            let inv=Option::<FieldElement>::from(inv_ct).expect("nonzero invert");

            let bytes=inv.to_bytes();

            let mut low=[0u8;8];

            low.copy_from_slice(&bytes[24..32]);

            checksum^=u64::from_be_bytes(low);
        }

        let dur=start.elapsed();

        println!(
            "k256_inv run {:02}: {:?}, checksum={}",
            run+1,dur,checksum
        );

        results.push(dur.as_secs_f64());
    }

    summarize("k256_inv summary_sec",&results);
}

fn main(){

    let runs=20usize;

    let iters=10_000usize;

    println!("== inverse only tournament ==");

    warmup();

    correctness_spot();

    println!("\n== standalone inverse ==");

    bench_inv_variant(InvChoice::Baseline,iters,runs);
    bench_inv_variant(InvChoice::Window3,iters,runs);
    bench_inv_variant(InvChoice::Window4,iters,runs);
    bench_inv_variant(InvChoice::Window5,iters,runs);
    bench_inv_variant(InvChoice::Window6,iters,runs);

    bench_inv_k256(iters,runs);

    println!("\n== inverse only tournament complete ==");
}
