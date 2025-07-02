extern crate f_psi;

use std::time::Instant;

use criterion::{criterion_group, criterion_main, Criterion};
use f_psi::psi::{self, Point, DIM, R, SIDE_LEN};
use f_psi::{protocol, psi_test};
use rand::Rng;

// 随机采集num个点
fn sample_test_data_points(num: usize) -> Vec<Point> {
    let mut rng = rand::thread_rng();
    let mut points: Vec<Point> = Vec::with_capacity(num);
    for _ in 0..num {
        let mut point: Point = [0u64; DIM];
        for i in 0..DIM {
            point[i] = rng.gen_range(SIDE_LEN..=(1 << 31));
        }
        points.push(point);
    }
    return points;
}

pub fn lp_single(b: &mut Criterion) {
    let n_r = 1;
    let n_s = 1;
    let metric = 1;

    let data_r = sample_test_data_points(n_r);
    let mut data_s = sample_test_data_points(n_s);

    data_s[0][0] = data_r[0][0] - R / 2;
    data_s[0][1] = data_r[0][1] + R / 2;

    println!("r: {:?}, s: {:?}", data_r, data_s);

    let mut psi_rec = psi_test::Receiver::new_single_lp(); // 创建Recv
    let psi_sed = psi_test::Sender::new_single_lp(1, psi_rec.publish_pk(), metric);

    b.bench_function("lp_single", |b| {
        b.iter(|| {
            let msg1 = psi_rec.lp_single_msg(&data_r, metric);

            let msg: (
                (
                    curve25519_dalek::RistrettoPoint,
                    curve25519_dalek::RistrettoPoint,
                ),
                std::collections::HashSet<u32>,
            ) = psi_sed.lp_send_msg_single(&msg1, &data_s[0], 0); // 发送单个消息

            let res = psi_rec.lp_single_post_process(&msg);
            println!("res: {:?}", res);
        })
    });
}

pub fn lin_single(b: &mut Criterion) {
    let n_r = 1;
    let n_s = 1;

    let data_r = sample_test_data_points(n_r);
    let mut data_s = sample_test_data_points(n_s);

    data_s[0][0] = data_r[0][0] - R / 2;
    data_s[0][1] = data_r[0][1] + R / 2;

    println!("r: {:?}, s: {:?}", data_r, data_s);

    let mut psi_rec = psi_test::Receiver::new_single_lin(); // 创建Recv
    let psi_sed = psi_test::Sender::new_single_lin(psi_rec.publish_pk());

    b.bench_function("lin_single", |b| {
        b.iter(|| {
            let msg1 = psi_rec.lin_single_msg(&data_r);

            let msg: Vec<(
                curve25519_dalek::RistrettoPoint,
                curve25519_dalek::RistrettoPoint,
            )> = psi_sed.lin_send_msg_single(&msg1, &data_s[0], 0); // 发送单个消息

            let res = psi_rec.lin_single_post_process(&msg);
            println!("res: {:?}", res);
        })
    });
}

pub fn lp_low(n_r: usize, n_s: usize) {
    println!("n: {}, m: {}, d:{}, delta:{}", n_r, n_s, DIM, R);
    let data_r = sample_test_data_points(n_r);
    let mut data_s = sample_test_data_points(n_s);
    data_s[9][0] = data_r[7][0] - R / 2;
    data_s[9][1] = data_r[7][1] + R / 2;
    data_s[11][0] = data_r[7][0] + R;
    data_s[11][1] = data_r[7][1] - R;
    let (psi_rcr, psi_sdr) = protocol::setup(n_r, n_s, true, 0);

    let now = Instant::now();
    protocol::run_standard_apart(psi_rcr, psi_sdr, data_r, data_s);
    let end = now.elapsed();
    println!("lin_low time: {:?}", end.as_millis());
}

pub fn lin_low(n_r: usize, n_s: usize) {
    println!("n: {}, m: {}, d:{}, delta:{}", n_r, n_s, DIM, R);
    let data_r = sample_test_data_points(n_r);
    let mut data_s = sample_test_data_points(n_s);
    data_s[9][0] = data_r[7][0] - R / 2;
    data_s[9][1] = data_r[7][1] + R / 2;
    data_s[11][0] = data_r[7][0] + R;
    data_s[11][1] = data_r[7][1] - R;

    let (psi_rcr, psi_sdr) = protocol::setup(n_r, n_s, true, 2);

    let now = Instant::now();
    protocol::run_standard_lp(psi_rcr, psi_sdr, data_r, data_s, 2);

    let end = now.elapsed();
    println!("lin_low time: {:?}", end.as_millis());
}

criterion_group!(
    name=benches;
    config = Criterion::default().significance_level(0.01).sample_size(10);
    targets=
    lp_single,
    lin_single,);

criterion_main!(benches);
