use crate::{okvs, psi}; // 引入 psi 模块

use std::io::{Read, Write};
use std::net::{TcpListener, TcpStream};
use std::sync::mpsc as std_mpsc; // 引入标准多生产者-单消费者通道
use std::thread;
use std::time::Instant; // 引入线程模块

// 设置函数，返回 psi::Receiver 和 psi::Sender
pub fn setup(
    dim: usize,
    radius: usize,
    num_n: usize,
    num_m: usize,
    apart: bool,
    sigma: bool,
) -> (psi::Receiver, psi::Sender) {
    let psi_rec = psi::Receiver::new(dim, radius as u64, num_m as u64, apart, sigma); // 创建Recv
    let psi_sed = psi::Sender::new(
        dim,
        radius as u64,
        num_n as u64,
        psi_rec.publish_pk(),
        apart,
        sigma,
    ); // 创建Sender

    return (psi_rec, psi_sed); // 返回Recv和Sender
}

// 标准 LP 运行函数
pub fn run_standard(
    mut psi_rec: psi::Receiver, // 可变Recv
    psi_sed: psi::Sender,       // Sender
    data_r: Vec<psi::Point>,    // 数据集 R
    data_s: Vec<psi::Point>,    // 数据集 S
    port: usize,                // 端口号
) -> (usize, usize, u128) {
    let (done_tx, done_rx) = std_mpsc::channel::<()>(); // 创建完成信号通道
    let (len_tx, len_rx) = std_mpsc::channel(); // 创建长度通道

    let now = Instant::now();

    let msg1 = psi_rec.inf_msg_apart(&data_r); // 获取 Lin 消息
    let msg1_com = bincode::serialize(&msg1).unwrap();
    let msg1_com_len = msg1_com.len();

    thread::spawn(move || {
        let mut stream =
            TcpStream::connect(format!("127.0.0.1:{}", port)).expect("Failed to connect to sender");

        // 发送 okvs 消息
        stream.write_all(&msg1_com).expect("Failed to send okvs"); // 发送 Lin 消息

        let mut count: u32 = 0u32; // 初始化计数器
        let mut buffer = Vec::new(); // 初始化缓冲区

        loop {
            let mut temp_buf = [0; 1024 * 1024];
            let n = stream.read(&mut temp_buf).unwrap();
            if n == 0 {
                break;
            } // 对端关闭
            buffer.extend_from_slice(&temp_buf[..n]);

            loop {
                if buffer.len() < 4 {
                    break;
                } // 长度前缀不足
                let msg_len = u32::from_le_bytes(buffer[..4].try_into().unwrap()) as usize;
                if buffer.len() < 4 + msg_len {
                    break;
                } // 消息不完整

                let single_msg_bytes = buffer[4..4 + msg_len].to_vec();
                let single_msg =
                    bincode::deserialize::<okvs::PointPair>(&single_msg_bytes).unwrap();

                // 处理消息
                count += psi_rec.inf_post_process_apart(&single_msg);

                // 移除已处理字节
                buffer.drain(..4 + msg_len);
            }

            // println!("count: {}, finished: {}", count, (count == 2)); // 打印接收消息数量
        }

        done_tx.send(()).expect("Failed to send done signal"); // 发送完成信号
    });

    thread::spawn(move || {
        // 绑定到指定端口
        let listener =
            TcpListener::bind(format!("127.0.0.1:{}", port)).expect("Failed to bind to address");

        let (mut stream, _) = listener.accept().expect("Failed to accept connection");

        // 接收并解码 okvs 消息
        let mut okvs_buf = vec![0u8; msg1_com_len];
        stream.read_exact(&mut okvs_buf).unwrap(); // 从流中读取数据
        let msg1_okvs: Vec<okvs::Encoding> = bincode::deserialize(&okvs_buf).unwrap(); // 反序列化消息

        let mut total_len = 0usize;
        for i in 0..data_s.len() {
            // 遍历数据集 S
            let single_msg = psi_sed.inf_send_msg_single_apart(&msg1_okvs, &data_s[i], i); // 发送单个 Lin 消息
            let single_msg_serialize = bincode::serialize(&single_msg).unwrap();
            let single_msg_len = single_msg_serialize.len();
            total_len += single_msg_len; // 累加消息长度

            stream
                .write_all(&(single_msg_len as u32).to_le_bytes())
                .unwrap(); // 发送长度
            stream.write_all(&single_msg_serialize).unwrap(); // 发送数据
        }

        len_tx.send(total_len).expect("Failed to send length"); // 发送总长度
    });

    // 等待任务完成
    done_rx.recv().expect("Failed to receive done signal"); // 接收完成信号

    let end = now.elapsed();

    let msg_rec_len = len_rx.recv().expect("Failed to receive length");

    (msg1_com_len, msg_rec_len, end.as_millis())
}
