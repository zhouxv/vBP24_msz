use core::time;
use std::time::Instant;
use std::{env, u128};

use f_psi::psi::{Point, DIM, R, SIDE_LEN};
use f_psi::{protocol, psi_test};
use rand::Rng;

fn main() {
    let args: Vec<String> = env::args().collect();
    println!("程序名称+参数: {:?}", args);

    let pro_type = args[1].parse::<usize>().unwrap_or(0);
    let times = args[2].parse::<usize>().unwrap_or(5);
    let pt_nums = args[3].parse::<usize>().unwrap_or(256);
    let metric = args[4].parse::<usize>().unwrap_or(0);

    match pro_type {
        0 => {
            // 低维协议
            test_all_low(times, pt_nums, metric);
        }
        1 => {
            // 1 vs 1
            test_all_single(times);
        }
        _ => {}
    }
}

pub fn test_all_single(times: usize) {
    println!("D:{}, R:{}", DIM, R);
    println!("-----------------无穷范数-------------------");

    let mut a = 0;
    let mut b = 0;
    let mut c = 0;
    for _ in 0..times {
        let (d, e, f) = fuzzy_macthing_infinity(false);
        println!("{} {} {}", d, e, f);
        a = a + d;
        b = b + e;
        c = c + f;
    }

    println!(
        "msg {} msg2 {} 总时间 {} μs",
        a / times,
        b / times,
        c / times as u128
    );

    println!("------------------L1---------------------");
    a = 0;
    b = 0;
    c = 0;
    for _ in 0..times {
        let (d, e, f) = fuzzy_macthing_lp(1, false);
        println!("{} {} {}", d, e, f);
        a = a + d;
        b = b + e;
        c = c + f;
    }

    println!(
        "msg {} msg2 {} 总时间 {} μs",
        a / times,
        b / times,
        c / times as u128
    );

    println!("------------------L2---------------------");

    a = 0;
    b = 0;
    c = 0;
    for _ in 0..times {
        let (d, e, f) = fuzzy_macthing_lp(2, false);
        println!("{} {} {}", d, e, f);
        a = a + d;
        b = b + e;
        c = c + f;
    }

    println!(
        "msg {} msg2 {} 总时间 {} μs",
        a / times,
        b / times,
        c / times as u128
    );
}

pub fn fuzzy_macthing_infinity(setting: bool) -> (usize, usize, u128) {
    let data_r = sample_test_data_points(1);
    let mut data_s = sample_test_data_points(1);
    if setting {
        for i in 0..DIM {
            data_s[0][i] = data_r[0][i] - R / 2;
        }
    }

    // println!("r: {:?}, s: {:?}", data_r, data_s);

    let mut psi_rec = psi_test::Receiver::new_single_lin(); // 创建Recv
    let psi_sed = psi_test::Sender::new_single_lin(psi_rec.publish_pk());

    let now = Instant::now();

    let msg1 = psi_rec.lin_single_msg(&data_r);
    let msg1_com = bincode::serialize(&msg1).unwrap();

    let msg: Vec<(
        curve25519_dalek::RistrettoPoint,
        curve25519_dalek::RistrettoPoint,
    )> = psi_sed.lin_send_msg_single(&msg1, &data_s[0], 0); // 发送单个消息

    //
    let msg_com = bincode::serialize(&msg).unwrap();

    let res = psi_rec.lin_single_post_process(&msg);

    //
    let end = now.elapsed();

    print!("res: {:?} ", res);

    (msg1_com.len(), msg_com.len(), end.as_millis())
}

pub fn fuzzy_macthing_lp(metric: usize, setting: bool) -> (usize, usize, u128) {
    let data_r = sample_test_data_points(1);
    let mut data_s = sample_test_data_points(1);

    if setting {
        for i in 0..DIM {
            data_s[0][i] = data_r[0][i] - 1;
        }
    }

    // println!("r: {:?}, s: {:?}", data_r, data_s);

    let mut psi_rec = psi_test::Receiver::new_single_lp(); // 创建Recv
    let psi_sed = psi_test::Sender::new_single_lp(1, psi_rec.publish_pk(), metric);

    let now = Instant::now();
    let msg1 = psi_rec.lp_single_msg(&data_r, metric);
    let msg1_com = bincode::serialize(&msg1).unwrap();

    let msg = psi_sed.lp_send_msg_single(&msg1, &data_s[0], 0); // 发送单个消息
    let msg_com = bincode::serialize(&msg).unwrap();

    let res = psi_rec.lp_single_post_process(&msg);
    let end = now.elapsed();

    print!("res: {:?} ", res);

    (msg1_com.len(), msg_com.len(), end.as_micros())
}

pub fn test_all_low(times: usize, pt_nums: usize, metric: usize) {
    if metric == 0 {
        lin_low(pt_nums, pt_nums, times);
    }

    if metric == 1 {
        lp_low(pt_nums, pt_nums, 1, times);
    }

    if metric == 2 {
        lp_low(pt_nums, pt_nums, 2, times);
    }
}

fn lin_low(n_r: usize, n_s: usize, times: usize) {
    println!("n: {}, m: {}, d:{}, R:{}", n_r, n_s, DIM, R);

    let mut a = 0;
    let mut b = 0;
    let mut c = 0;

    for _ in 0..times {
        let data_r = sample_test_data_points(n_r);
        let mut data_s = sample_test_data_points(n_s);

        for i in 0..DIM {
            data_s[9][i] = data_r[7][i] - R / 2;
            data_s[11][i] = data_r[7][i] + R;
        }

        let (psi_rcr, psi_sdr) = protocol::setup(n_r, n_s, true, 0);

        let (d, e, f) = protocol::run_standard_apart(psi_rcr, psi_sdr, data_r, data_s);

        println!("msg1 {} msg2 {} times {}", d, e, f);
        a = a + d;
        b = b + e;
        c = c + f;
    }

    let avg_a = a / times;
    let avg_b = b / times;
    let avg_c = c / times as u128;
    let avg_com = (avg_a + avg_b) as f64 / 1024.0 / 1024.0;

    println!(
        "{} 字节 {} 字节 {} 毫秒 {} MB",
        avg_a, avg_b, avg_c, avg_com
    );
}

fn lp_low(n_r: usize, n_s: usize, metric: usize, times: usize) {
    println!(
        "n: {}, m: {}, D:{}, R:{}, metric:{}",
        n_r, n_s, DIM, R, metric
    );

    let mut a = 0;
    let mut b = 0;
    let mut c = 0;

    for _ in 0..times {
        let data_r = sample_test_data_points(n_r);
        let mut data_s = sample_test_data_points(n_s);

        for i in 0..DIM {
            data_s[9][i] = data_r[7][i] - R / 2;
            data_s[11][i] = data_r[7][i] + R;
        }

        let (psi_rcr, psi_sdr) = protocol::setup(n_r, n_s, true, metric);

        let (d, e, f) = protocol::run_standard_lp(psi_rcr, psi_sdr, data_r, data_s, metric as u32);

        println!("msg1 {} msg2 {} times {}", d, e, f);

        a = a + d;
        b = b + e;
        c = c + f;
    }

    let avg_a = a / times;
    let avg_b = b / times;
    let avg_c = c / times as u128;
    let avg_com = (avg_a + avg_b) as f64 / 1024.0 / 1024.0;

    println!(
        "{} 字节 {} 字节 {} 毫秒 {} MB",
        avg_a, avg_b, avg_c, avg_com
    );
}

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
