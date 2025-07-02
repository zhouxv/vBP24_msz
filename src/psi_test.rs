use std::collections::HashSet;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

use crate::okvs;
use crate::psi::Point;
use crate::psi::DIM;
use crate::psi::R;
use crate::psi::R_L2;
use fxhash::hash64;

pub struct Receiver {
    window: usize,                         // 窗口大小
    okvs_size: usize,                      // okvs 大小
    pk: RistrettoPoint,                    // 公钥
    sk: Scalar,                            // 私钥
    _pre_data: Vec<Vec<(Scalar, Scalar)>>, // 预处理数据
    _okvsgen: Vec<okvs::OkvsGen>,          // okvs生成器
}

impl Receiver {
    // 返回公钥
    pub fn publish_pk(&self) -> RistrettoPoint {
        return self.pk;
    }

    // 返回窗口大小
    pub fn get_windowsize(&self) -> usize {
        return self.window;
    }

    // 1 vs 1 test new L∞
    pub fn new_single_lin() -> Self {
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
        let sk: Scalar = Scalar::random(&mut rng); // 随机生成私钥
        let pk: RistrettoPoint = &sk * RISTRETTO_BASEPOINT_TABLE; // 根据私钥计算公钥
        let mut _pre_data: Vec<Vec<(Scalar, Scalar)>> = Vec::with_capacity(DIM); // 初始化预处理数据
        let mut _okvsgen: Vec<okvs::OkvsGen> = Vec::with_capacity(DIM); // 初始化 okvs 生成器

        let window = 2 * R + 1; // 否则窗口大小为 (2 * R + 1)

        for _ in 0..DIM {
            let mut pair: Vec<(Scalar, Scalar)> = Vec::with_capacity(window as usize); // 为每一维度初始化一个空的预处理数据
            for _ in 0..window {
                let tem: Scalar = Scalar::random(&mut rng); // 随机生成数据
                pair.push((tem, tem * sk)); // 将数据和数据乘以私钥存入配对中
            }
            _pre_data.push(pair); // 将每一维度的预处理数据加入 _pre_data
            _okvsgen.push(okvs::OkvsGen::new(window)); // 为每一维度初始化 okvs 生成器
        }

        return Receiver {
            window: window as usize,
            okvs_size: window as usize,
            pk,
            sk,
            _pre_data,
            _okvsgen,
        };
    }

    /*
    为每个维度处理一组点集，并将结果编码为 okvs::Encoding 类型的向量
    */
    pub fn lin_single_msg(&mut self, pt_set: &Vec<Point>) -> Vec<okvs::Encoding> {
        let mut result: Vec<okvs::Encoding> = Vec::with_capacity(DIM); // 初始化OKVS结果向量
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new(); // 初始化列表，存储键值对
                                                                 // 遍历每一个维度
        for i in 0..DIM {
            // 遍历每个接收者的点 pt
            // 将每个 Point 与 对应的预处理数目 对应
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            {
                let min: u64 = pt[i] - R as u64;
                for (j, pre_val) in pre_window.iter().enumerate() {
                    let key_ij = hash64(&(min + j as u64)); // 计算每个点的哈希值
                    list.push((key_ij, *pre_val)); // 将键值对添加到列表中
                }
            }
            result.push(self._okvsgen[i].encode(&list)); // 对每个维度的数据进行编码
            list.clear(); // 清空列表
        }
        self._okvsgen.clear(); // 清空 okvs 生成器
        self._pre_data.clear(); // 清空预处理数据

        // println!("{} {}", result.len(), result[0].len());

        return result; // 返回编码结果
    }

    // 后处理操作，检查消息是否有效
    #[inline]
    pub fn lin_single_post_process(&mut self, msg_sender: &okvs::Encoding) -> u32 {
        for (u, v) in msg_sender.iter() {
            if (self.sk * u) == *v {
                return 1; // 如果匹配，返回 1
            }
        }
        return 0; // 否则返回 0
    }

    pub fn new_single_lp() -> Self {
        let num_item: u64 = 1;
        // 生成随机的私钥和公钥
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
        let sk: Scalar = Scalar::random(&mut rng); // 随机生成私钥
        let pk: RistrettoPoint = &sk * RISTRETTO_BASEPOINT_TABLE; // 根据私钥计算公钥
        let mut _pre_data: Vec<Vec<(Scalar, Scalar)>> = Vec::with_capacity(DIM); // 初始化预处理数据
        let mut _okvsgen: Vec<okvs::OkvsGen> = Vec::with_capacity(DIM); // 初始化 okvs 生成器

        let window: usize = (2 * R + 1) as usize;

        let n: u64 = num_item * window as u64; // 计算输出的大小
        for _ in 0..DIM {
            let mut pair: Vec<(Scalar, Scalar)> = Vec::with_capacity(n as usize); // 为每一维度初始化一个空的预处理数据
            for _ in 0..n {
                let tem: Scalar = Scalar::random(&mut rng); // 随机生成数据
                pair.push((tem, tem * sk)); // 将数据和数据乘以私钥存入配对中
            }
            _pre_data.push(pair); // 将每一维度的预处理数据加入 _pre_data
            _okvsgen.push(okvs::OkvsGen::new(n)); // 为每一维度初始化 okvs 生成器
        }

        return Receiver {
            window,
            okvs_size: n as usize,
            pk,
            sk,
            _pre_data,
            _okvsgen,
        };
    }

    pub fn lp_single_msg(&mut self, pt_set: &Vec<Point>, metric: usize) -> Vec<okvs::Encoding> {
        let mut result: Vec<okvs::Encoding> = Vec::with_capacity(DIM); // 初始化结果向量
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new(); // 初始化列表，存储键值对
        let possible_vals = 2 * R as usize + 1; // 计算可能的值
                                                // 遍历每一个维度
        for i in 0..DIM {
            // 遍历每个接收者的点 pt
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            // 使用滑动窗口
            {
                let min = pt[i] - R as u64;
                for j in 0..possible_vals {
                    let key_ij = hash64(&(min + j as u64)); // 计算每个点的哈希值
                    let tem: (Scalar, Scalar) = pre_window[j];
                    let mut diff_abs = if j as u64 > R {
                        j as u64 - R
                    } else {
                        R - j as u64
                    };
                    if metric == 2 {
                        diff_abs *= diff_abs; // 如果是 L2 距离，计算平方
                    }
                    // tem.0=a, tem.1=a*sk
                    // diff_abs=|x-R|^p
                    list.push((key_ij, (tem.0, tem.1 + Scalar::from(diff_abs))));
                    // 添加到列表
                }
            }
            result.push(self._okvsgen[i].encode(&list)); // 对每个维度的数据进行编码
            list.clear(); // 清空列表
        }
        self._okvsgen.clear(); // 清空 okvs 生成器
        self._pre_data.clear(); // 清空预处理数据
        return result; // 返回编码结果
    }

    // 后处理操作，检查消息是否有效（适用于 LP 模式）
    #[inline]
    pub fn lp_single_post_process(&mut self, msg_sender: &(okvs::PointPair, HashSet<u32>)) -> u32 {
        let x = hash64(
            &(msg_sender.0 .1 - self.sk * msg_sender.0 .0)
                .compress()
                .to_bytes(),
        ) as u32; // 计算哈希值

        // println!("x2: {}", x);

        if msg_sender.1.contains(&x) {
            return 1; // 如果哈希值匹配，返回 1
        }
        return 0; // 否则返回 0
    }
}

pub struct Sender {
    m: u64,                                                // 发送者的消息数量
    window: usize,                                         // 窗口大小
    pk: RistrettoPoint,                                    // 公钥
    _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)>, // 存储coin的向量，包含 RistrettoPoint 和 Scalar
    _coins_lp_set: Vec<HashSet<u32>>,                      // 存储 Lp 模式下的coin集合
    _coins_lp: Vec<RistrettoPoint>,                        // 存储 Lp 模式下的coin
}

impl Sender {
    // 1 vs 1 test new
    pub fn new_single_lin(pk: RistrettoPoint) -> Self {
        let m: u64 = 1; // 计算消息数量
        let mut rng = rand::thread_rng(); // 创建随机数生成器
        let mut _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)> =
            Vec::with_capacity(m as usize); // 初始化coin向量

        for _ in 0..m {
            // 生成coin
            let a = Scalar::random(&mut rng); // 随机生成一个 Scalar
            let b = Scalar::random(&mut rng); // 随机生成另一个 Scalar

            // _coins.0=a,_coins.1=a*pk,_coins.2=b
            _coins.push((&a * RISTRETTO_BASEPOINT_TABLE, a * pk, b)); // 将coin添加到向量中
        }

        return Sender {
            // 返回 Sender 实例
            m,
            window: 1,
            pk,
            _coins,
            _coins_lp_set: Vec::new(),
            _coins_lp: Vec::new(),
        };
    }

    // 发送单个消息
    #[inline]
    pub fn lin_send_msg_single(
        &self,
        encodings: &Vec<okvs::Encoding>, // 编码向量
        pt: &Point,                      // 点
        index: usize,                    // 索引
    ) -> okvs::Encoding {
        let coin_window = &self._coins[index * self.window..(index + 1) * self.window]; // 获取当前窗口的coin
        let mut uv: okvs::Encoding =
            vec![(RistrettoPoint::identity(), RistrettoPoint::identity()); 1]; // 初始化结果编码
        let mut tem: okvs::PointPair; // 临时变量

        // 遍历每个可能的块
        // 因为在recv计算用的是对应Point的block, sender则不知道，到底是哪个block
        for (i, coins) in coin_window.iter().enumerate() {
            for j in 0..DIM {
                let key_ij = hash64(&(pt[j] as u64)); // 计算每个点的哈希值
                tem = okvs::okvs_decode(&encodings[j], key_ij); // 解码
                uv[i].0 += tem.0; // 更新结果
                uv[i].1 += tem.1; // 更新结果
            }
            // 完成计算
            // uv[i].0= b*x+a
            // uv[i].1= b*(x*sk)+a*pk
            uv[i].0 = coins.2 * uv[i].0 + coins.0; // 最终计算
            uv[i].1 = coins.2 * uv[i].1 + coins.1; // 最终计算
        }
        return uv; // 返回编码结果
    }

    // 构造函数，初始化 Sender 实例
    pub fn new_single_lp(num_item: u64, pk: RistrettoPoint, metric: usize) -> Self {
        let window: usize = 1; // 根据 apart 参数设置窗口大小
        let m = num_item * window as u64; // 计算消息数量
        let mut rng = rand::thread_rng(); // 创建随机数生成器
        let mut _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)> =
            Vec::with_capacity(m as usize); // 初始化coin向量

        let metric_window = if metric == 2 {
            // 根据度量选择窗口大小
            R_L2 as usize + 1
        } else {
            R as usize + 1
        };

        let mut _coins_lp_set: Vec<HashSet<u32>> = Vec::with_capacity(m as usize); // 初始化 Lp 模式下的coin集合
        let mut _coins_lp: Vec<RistrettoPoint> = Vec::with_capacity(m as usize); // 初始化 Lp 模式下的coin
        for _ in 0..m {
            // 生成coin
            let a = Scalar::random(&mut rng); // 随机生成一个 Scalar
            let b = Scalar::random(&mut rng); // 随机生成另一个 Scalar
            let c = Scalar::random(&mut rng); // 随机生成第三个 Scalar

            // _coins.0=a,_coins.1=a*pk,_coins.2=b
            _coins.push((&a * RISTRETTO_BASEPOINT_TABLE, a * pk, b)); // 将coin添加到向量中
            _coins_lp.push(&c * RISTRETTO_BASEPOINT_TABLE); // 将 Lp coin添加到向量中
            let mut hashtab: HashSet<u32> = HashSet::with_capacity(metric_window); // 初始化哈希集合
            for i in 0..metric_window as u64 {
                // 计算哈希值并存储
                // g_ j= c+b*i
                let g_j = &(c + b * Scalar::from(i)) * RISTRETTO_BASEPOINT_TABLE; // 计算 g_j
                let x = hash64(&g_j.compress().to_bytes()) as u32;
                // println!("x: {}", x);
                hashtab.insert(x); // 将哈希值插入集合
            }
            _coins_lp_set.push(hashtab); // 将哈希集合添加到 Lp coin集合中
        }

        return Sender {
            // 返回 Sender 实例
            m,
            window,
            pk,
            _coins,
            _coins_lp_set,
            _coins_lp,
        };
    }

    // 发送单个消息（适用于 LP 模式）
    #[inline]
    pub fn lp_send_msg_single(
        &self,
        encodings: &Vec<okvs::Encoding>, // 编码向量
        pt: &Point,                      // 点
        index: usize,                    // 索引
    ) -> (okvs::PointPair, HashSet<u32>) {
        let coins = &self._coins[index]; // 获取当前窗口的coin
        let mut uv: okvs::PointPair = (RistrettoPoint::identity(), RistrettoPoint::identity()); // 初始化结果
        let mut tem: okvs::PointPair; // 临时变量

        for j in 0..DIM {
            let key_ij = hash64(&pt[j]); // 计算每个点的哈希值
            tem = okvs::okvs_decode(&encodings[j], key_ij); // 解码
            uv.0 += tem.0; // 更新结果
            uv.1 += tem.1; // 更新结果
        }

        // 完成计算
        uv.0 = coins.2 * uv.0 + coins.0; // 最终计算
        uv.1 = coins.2 * uv.1 + coins.1 + &self._coins_lp[index]; // 最终计算

        return (uv, self._coins_lp_set[index].clone()); // 返回结果和 Lp coin集合
    }

    // 返回输出的大小
    pub fn get_output_size(&self) -> u64 {
        return self.m; // 返回消息数量
    }

    // 返回窗口大小
    pub fn get_windowsize(&self) -> usize {
        return self.window; // 返回窗口大小
    }
}
