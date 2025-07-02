use crate::psi; // 引入 psi 模块

use std::sync::mpsc as std_mpsc; // 引入标准多生产者-单消费者通道
use std::thread;
use std::time::Instant; // 引入线程模块

// 设置函数，返回 psi::Receiver 和 psi::Sender
pub fn setup(
    num_n: usize,
    num_m: usize,
    apart: bool,
    metric: usize,
) -> (psi::Receiver, psi::Sender) {
    let psi_rec = psi::Receiver::new(num_n as u64, apart); // 创建Recv
    let psi_sed = psi::Sender::new(num_m as u64, psi_rec.publish_pk(), apart, metric as u32); // 创建Sender
    return (psi_rec, psi_sed); // 返回Recv和Sender
}

// 标准apart运行函数
pub fn run_standard_apart(
    mut psi_rec: psi::Receiver, // Recv
    psi_sed: psi::Sender,       // Sender
    data_r: Vec<psi::Point>,    // 数据集 R
    data_s: Vec<psi::Point>,    // 数据集 S
) -> (usize, usize, u128) {
    let (done_tx, done_rx) = std_mpsc::channel::<()>(); // 创建完成信号通道
    let (sender, receiver) = std_mpsc::channel(); // 创建发送和接收通道
    let (len_tx, len_rx) = std_mpsc::channel(); // 创建长度通道

    let now = Instant::now();

    let msg1 = psi_rec.msg_apart(&data_r); // 获取apart的消息

    let msg1_com = bincode::serialize(&msg1).unwrap();
    let msg1_com_len = msg1_com.len();

    // 发送线程
    thread::spawn(move || {
        for i in 0..data_s.len() {
            // 遍历数据集 S
            let msg = psi_sed.send_msg_single_apart(&msg1, &data_s[i], i); //

            if i == 0 {
                let msg2_com = bincode::serialize(&msg).unwrap();
                let len = msg2_com.len() * data_s.len();

                len_tx.send(len).expect("Failed to send length");
            }

            if sender.send(msg).is_err() {
                // 如果发送失败
                println!("Receiver has been dropped!"); // 打印错误信息
                break; // 退出循环
            }
        }
    });

    // 接收线程
    thread::spawn(move || {
        let mut count = 0u32; // 初始化计数器
        for msg2 in receiver.iter() {
            // 遍历接收到的消息
            count += psi_rec.post_process_apart(&msg2); // 处理消息并更新计数
        }
        // println!("count: {}", count); // 打印计数
        done_tx.send(()).expect("Failed to send done signal"); // 发送完成信号
    });

    // 等待任务完成
    done_rx.recv().expect("Failed to receive done signal"); // 接收完成信号

    let end = now.elapsed();

    // -----------------------------------------------------------------------------------------------------------

    let msg_rec_len = len_rx.recv().expect("Failed to receive length");

    (msg1_com_len, msg_rec_len, end.as_millis())
}

// 标准 LP 运行函数
pub fn run_standard_lp(
    mut psi_rec: psi::Receiver, // 可变Recv
    psi_sed: psi::Sender,       // Sender
    data_r: Vec<psi::Point>,    // 数据集 R
    data_s: Vec<psi::Point>,    // 数据集 S
    metric: u32,                // 度量
) -> (usize, usize, u128) {
    let (done_tx, done_rx) = std_mpsc::channel::<()>(); // 创建完成信号通道
    let (sender, receiver) = std_mpsc::channel(); // 创建发送和接收通道
    let (len_tx, len_rx) = std_mpsc::channel(); // 创建长度通道

    let now = Instant::now();

    let msg1 = psi_rec.lp_msg_apart(&data_r, metric); // 获取 LP 消息

    let msg1_com = bincode::serialize(&msg1).unwrap();
    let msg1_com_len = msg1_com.len();

    // 发送线程
    thread::spawn(move || {
        for i in 0..data_s.len() {
            // 遍历数据集 S
            let msg = psi_sed.lp_send_msg_single_apart(&msg1, &data_s[i], i); // 发送单个 LP 消息

            if i == 0 {
                let msg2_com = bincode::serialize(&msg).unwrap();
                let len = msg2_com.len() * data_s.len();
                len_tx.send(len).expect("Failed to send length");
            }

            if sender.send(msg).is_err() {
                // 如果发送失败
                println!("Receiver has been dropped!"); // 打印错误信息
                break; // 退出循环
            }
        }
    });

    // 接收线程
    thread::spawn(move || {
        let mut count = 0u32; // 初始化计数器
        for msg2 in receiver.iter() {
            // 遍历接收到的消息
            count += psi_rec.lp_post_process_apart(&msg2); // 处理消息并更新计数
        }
        // println!("Lp metric {}, count: {}", metric, count); // 打印度量和计数
        done_tx.send(()).expect("Failed to send done signal"); // 发送完成信号
    });

    // 等待任务完成
    done_rx.recv().expect("Failed to receive done signal"); // 接收完成信号

    let end = now.elapsed();

    let msg_rec_len = len_rx.recv().expect("Failed to receive length");

    (msg1_com_len, msg_rec_len, end.as_millis())
}

#[cfg(test)]
mod tests {

    use std::time::Instant;

    use crate::{
        protocol,
        psi::{self, Point, DIM, R, SIDE_LEN},
    };
    use rand::Rng;

    #[test]
    pub fn test_all() {
        for _ in 0..10 {
            lin_low(16, 16);
            lin_low(256, 256);
            lin_low(4096, 4096);

            lp_low(16, 16, 1);
            lp_low(256, 256, 1);
            lp_low(4096, 4096, 1);

            lp_low(16, 16, 1);
            lp_low(256, 256, 2);
            lp_low(4096, 4096, 2);
        }
    }

    #[test]
    pub fn lin_low_16() {
        lin_low(16, 16);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lin_low_256() {
        lin_low(256, 256);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lin_low_4096() {
        lin_low(4096, 4096);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lp_low_16_1() {
        lp_low(16, 16, 1);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lp_low_256_1() {
        lp_low(256, 256, 1);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lp_low_4096_1() {
        lp_low(4096, 4096, 1);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lp_low_16_2() {
        lp_low(16, 16, 2);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lp_low_256_2() {
        lp_low(256, 256, 2);
        println!("---------------------------------------");
    }

    #[test]
    pub fn lp_low_4096_2() {
        lp_low(4096, 4096, 2);
        println!("---------------------------------------");
    }

    fn lp_low(n_r: usize, n_s: usize, metric: usize) {
        println!("n: {}, m: {}, d:{}, delta:{}", n_r, n_s, DIM, R);
        let data_r = sample_test_data_points(n_r);
        let mut data_s = sample_test_data_points(n_s);
        data_s[9][0] = data_r[7][0] - R / 2;
        data_s[9][1] = data_r[7][1] + R / 2;
        data_s[11][0] = data_r[7][0] + R;
        data_s[11][1] = data_r[7][1] - R;
        let (psi_rcr, psi_sdr) = protocol::setup(n_r, n_s, true, metric);

        let now = Instant::now();
        protocol::run_standard_apart(psi_rcr, psi_sdr, data_r, data_s);
        let end = now.elapsed();
        println!("lin_low time: {:?} μs", end.as_micros());
    }

    fn lin_low(n_r: usize, n_s: usize) {
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
        println!("lin_low time: {:?} μs", end.as_micros());
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
}
