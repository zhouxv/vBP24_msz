#![allow(dead_code)]

extern crate f_psi;

use curve25519_dalek::scalar::Scalar;
use f_psi::okvs;
use fxhash::hash64;

use criterion::{criterion_group, criterion_main, Criterion};

// OKVS Bench!
pub fn bench1(b: &mut Criterion) {
    let n = 327680;
    let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
    println!("{} items, OkvsGen.Decode", n);
    for j in 0..n {
        list.push((hash64(&j), (Scalar::ONE, Scalar::ONE)));
    }

    b.bench_function("bench_okvs_decode", |b| {
        b.iter(|| {
            let mut okvs_instance = okvs::OkvsGen::new(n);
            let data = okvs_instance.encode(&list);
            let keys: Vec<u64> = (0..n).collect();
            okvs::okvs_decode_batch(&data, &keys);
        })
    });
}

criterion_group!(
    name=benches;
    config = Criterion::default().significance_level(0.01).sample_size(10);
    targets=
    bench1,);
criterion_main!(benches);
