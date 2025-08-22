use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

use crate::okvs;
use fxhash::hash64;

pub type Point = Vec<u64>;

/*
计算一个给定点 p 在一个 d 维空间中根据边长 sidele 所对应的块的左下角坐标 Point
*/
#[inline]
fn cell(p: &Point, dim: usize, radius: u64, sigma: bool) -> Point {
    let mut bot_left_corner: Point = vec![0u64; dim];
    let sidele = if sigma { radius * 4 } else { radius * 2 }; // 如果 sigma 为 true，则边长乘以 2
    for i in 0..dim {
        bot_left_corner[i] = p[i] / sidele;
    }
    return bot_left_corner;
}

/*
计算一个给定点 p, 根据一个半径 radius 和边长 sidele 所对应的块的左下角坐标 Point
计算边界点所在的块
*/
#[inline]
fn block(p: &Point, dim: usize, radius: u64, sigma: bool) -> Point {
    let mut min: Point = vec![0u64; dim];
    for i in 0..dim {
        // 计算给定点 p 在某个维度 i 上的坐标减去一个半径值
        // 确定一个新的坐标，这个坐标通常用于定义一个区域或块的边界
        min[i] = p[i] - radius;
    }
    return cell(&min, dim, radius, sigma);
}

// 计算两个点 p1 和 p2 之间的 L2 距离（欧几里得距离的平方）
#[inline]
fn l2_dist(p1: &Point, p2: &Point, dim: usize) -> u64 {
    let mut sum: u64 = 0;
    let mut diff: u64;
    for i in 0..dim {
        // 计算每个维度的差值的绝对值
        diff = (p1[i] as i64 - p2[i] as i64).abs() as u64;
        // 将差值的平方累加到总和
        sum += diff * diff;
    }
    return sum;
}

// 计算每个维度的差值的绝对值，并累加
#[inline]
fn l1_dist(p1: &Point, p2: &Point, dim: usize) -> u64 {
    let mut sum: u64 = 0;
    for i in 0..dim {
        sum += (p1[i] as i64 - p2[i] as i64).abs() as u64;
    }
    return sum;
}

// 计算两个点之间的无穷范数（切比雪夫距离）
#[inline]
fn l_inf_dist(p1: &Point, p2: &Point, dim: usize) -> u64 {
    let mut max_diff: u64 = 0; // 初始化最大差值
    for i in 0..dim {
        let diff = (p1[i] as i64 - p2[i] as i64).abs() as u64; // 计算每个维度的差值
        if diff > max_diff {
            max_diff = diff; // 更新最大差值
        }
    }
    return max_diff; // 返回无穷范数
}

// 计算给定点 p 相对于源点 source 的位置索引
#[inline]
fn get_position(p: &Point, source: &Point, dim: usize) -> usize {
    let mut pos: usize = 0;
    for i in 0..dim {
        // 如果 p 的某一维度大于 source 的对应维度，则更新位置索引
        if p[i] > source[i] {
            pos += 1 << i;
        }
    }
    return pos;
}

#[inline]
fn intersection(recv: &Receiver, p: &Point, metric: u32) -> Vec<Point> {
    let dim = recv.dim;
    let radius = recv.radius;
    let blk_cells = recv.blk_cells;
    let sigma = recv.sigma;

    let side_len = if sigma {
        radius * 4 // 如果 sigma 为 true，则边长乘以 4
    } else {
        radius * 2 // 否则，边长乘以 2
    };

    // 初始化结果向量，容量为 blk_cells
    let mut result: Vec<Point> = Vec::with_capacity(blk_cells);
    // 计算给定点 p 所在块的左下角坐标
    let blk = block(p, dim, radius, sigma);
    // 初始化交叉点
    let mut cross_point: Point = vec![0u64; dim];

    // 计算交叉点的坐标, 交叉点是2*sigma的单元格的右上角的点
    // cross_point是右上角
    for j in 0..dim {
        cross_point[j] = (blk[j] + 1) * side_len;
    }

    // let dist = match metric {
    //     2 => l2_dist(p, &cross_point, dim),
    //     1 => l1_dist(p, &cross_point, dim),
    //     0 => l_inf_dist(p, &cross_point, dim),
    //     _ => panic!("metric should be L1 or L2"),
    // };

    // 获取 cross_point 相对于源点 p 的位置索引
    let pos_ind = get_position(&cross_point, p, dim);

    // 遍历所有块
    for i in 0..blk_cells {
        let mut tem: Point = vec![0u64; dim];

        // // 根据度量选择半径
        // let r_lp = match metric {
        //     2 => recv.r_l2,
        //     1 => radius,
        //     0 => radius,
        //     _ => panic!("metric should be 0, 1, or 2"),
        // };

        // 当查询点 p 靠近块的左下角时，右上角的邻接块可能不在邻域内，可以直接跳过
        // if (dist > r_lp * 2) && (i == pos_ind) {
        //     continue;
        // }

        // 计算当前块的坐标
        for j in 0..dim {
            if (i >> j) & 1 == 1 {
                tem[j] = blk[j] + 1;
            } else {
                tem[j] = blk[j];
            }
        }
        // 将当前块的坐标添加到结果中
        result.push(tem);
    }
    // 返回结果
    return result;
}

pub struct Receiver {
    dim: usize,                            // 维度
    radius: u64,                           // 半径
    sigma: bool,                           // 是否使用 sigma
    window: usize,                         // 窗口大小
    n: u64,                                // 输出的大小
    pk: RistrettoPoint,                    // 公钥
    sk: Scalar,                            // 私钥
    _pre_data: Vec<Vec<(Scalar, Scalar)>>, // 预处理数据
    _okvsgen: Vec<okvs::OkvsGen>,          // okvs生成器
    r_l2: u64,
    side_len: u64,    // 边长，直径
    blk_cells: usize, // 块单元格数量
}

impl Receiver {
    // 构造函数，初始化 Receiver 实例
    pub fn new(dim: usize, radius: u64, num_item: u64, apart: bool, sigma: bool) -> Self {
        let blk_cells: usize = 1 << dim;
        let side_len = 2 * radius;
        let r_l2 = radius * radius;

        // 生成随机的私钥和公钥
        let mut rng: rand::rngs::ThreadRng = rand::thread_rng();
        let sk: Scalar = Scalar::random(&mut rng); // 随机生成私钥
        let pk: RistrettoPoint = &sk * RISTRETTO_BASEPOINT_TABLE; // 根据私钥计算公钥
        let mut _pre_data: Vec<Vec<(Scalar, Scalar)>> = Vec::with_capacity(dim); // 初始化预处理数据
        let mut _okvsgen: Vec<okvs::OkvsGen> = Vec::with_capacity(dim); // 初始化 okvs 生成器

        let window: usize = if apart == true {
            blk_cells * (2 * radius + 1) as usize // 如果 apart 为 true，窗口大小为 blk_cells * (2 * radius + 1)
        } else {
            (2 * radius + 1) as usize // 否则窗口大小为 (2 * radius + 1)
        };

        let n: u64 = num_item * window as u64; // 计算输出的大小

        let tem: Scalar = Scalar::random(&mut rng); // 随机生成数据
        let tem_sk = tem * sk; // 将数据和数据乘以私钥存入配对中

        let mut pair: Vec<(Scalar, Scalar)> = Vec::with_capacity(n as usize); // 为每一维度初始化一个空的预处理数据

        for _ in 0..n {
            // let tem: Scalar = Scalar::random(&mut rng); // 随机生成数据
            // pair.push((tem, tem * sk)); // 将数据和数据乘以私钥存入配对中
            pair.push((tem, tem_sk));
        }

        for _ in 0..dim {
            // let mut pair: Vec<(Scalar, Scalar)> = Vec::with_capacity(n as usize); // 为每一维度初始化一个空的预处理数据

            // for _ in 0..n {
            //     // let tem: Scalar = Scalar::random(&mut rng); // 随机生成数据
            //     // pair.push((tem, tem * sk)); // 将数据和数据乘以私钥存入配对中
            //     pair.push((tem, tem_sk));
            // }
            _pre_data.push(pair.clone()); // 将每一维度的预处理数据加入 _pre_data

            _okvsgen.push(okvs::OkvsGen::new(n)); // 为每一维度初始化 okvs 生成器
        }

        // println!("recv init end");

        return Receiver {
            dim,
            radius,
            sigma,
            window,
            n,
            pk,
            sk,
            _pre_data,
            _okvsgen,
            r_l2,
            side_len,
            blk_cells,
        };
    }

    pub fn inf_msg_apart(&mut self, pt_set: &Vec<Point>) -> Vec<okvs::Encoding> {
        let mut result: Vec<okvs::Encoding> = Vec::with_capacity(self.dim); // 初始化结果向量
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new(); // 初始化列表，存储键值对
        let possible_vals = 2 * self.radius as usize + 1; // 计算可能的值
                                                          // 遍历每一个维度
        for i in 0..self.dim {
            // 遍历每个接收者的点 pt
            for (pt, pre_window) in pt_set
                .iter()
                .zip(self._pre_data[i].windows(self.window).step_by(self.window))
            // 使用滑动窗口
            {
                // 遍历 [2R+1] * blk_cells 范围内的每个值
                let cels = intersection(self, pt, 0); // 获取交叉点
                                                      // blk_cells*possible_vals*n
                for k in 0..self.blk_cells {
                    let key;
                    if k >= cels.len() {
                        key = rand::random::<u64>(); // 随机生成 key
                    } else {
                        key = hash64(&cels[k]); // 根据交叉点计算 key
                    }
                    let min = pt[i] - self.radius as u64;

                    for j in 0..possible_vals {
                        let key_ij = hash64(&(min + j as u64)); // 计算每个点的哈希值
                        list.push((key ^ key_ij, pre_window[k * possible_vals + j]));
                        // 添加到列表
                    }
                }
            }
            result.push(self._okvsgen[i].encode(&list)); // 对每个维度的数据进行编码
            list.clear(); // 清空列表
        }
        self._okvsgen.clear(); // 清空 okvs 生成器
        self._pre_data.clear(); // 清空预处理数据
        return result; // 返回编码结果
    }

    // 后处理操作，检查消息是否有效（适用于 apart 模式）
    #[inline]
    pub fn inf_post_process_apart(&mut self, msg_sender: &okvs::PointPair) -> u32 {
        if (self.sk * msg_sender.0) == msg_sender.1 {
            return 1; // 如果匹配，返回 1
        }
        return 0; // 否则返回 0
    }

    // 返回公钥
    pub fn publish_pk(&self) -> RistrettoPoint {
        return self.pk;
    }

    // 刷新 Receiver 实例
    pub fn refresh(&mut self) {
        let mut rng = rand::thread_rng(); // 创建随机数生成器
        self._okvsgen.clear(); // 清空 okvs 生成器
        self._pre_data.clear(); // 清空预处理数据
        for _ in 0..self.dim {
            let mut pair = Vec::with_capacity(self.n as usize); // 为每一维度重新初始化预处理数据
            for _ in 0..self.n {
                let tem = Scalar::random(&mut rng); // 随机生成数据
                pair.push((tem, tem * self.sk)); // 存储数据和数据乘以私钥
            }
            self._pre_data.push(pair); // 更新预处理数据
            self._okvsgen.push(okvs::OkvsGen::new(self.n)); // 更新 okvs 生成器
        }
    }

    pub fn print(&self) {
        println!(
            "Receiver: dim={}, radius={}, sigma={}, window={}, n={}",
            self.dim, self.radius, self.sigma, self.window, self.n
        );
        println!(
            "r_l2: {}, side_len: {}, blk_cells: {}",
            self.r_l2, self.side_len, self.blk_cells
        );
    }
}

pub struct Sender {
    dim: usize,  // 维度
    radius: u64, // 半径
    sigma: bool,
    m: u64,                                                // 发送者的消息数量
    window: usize,                                         // 窗口大小
    pk: RistrettoPoint,                                    // 公钥
    _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)>, // 存储coin的向量，包含 RistrettoPoint 和 Scalar
    r_l2: u64,
    side_len: u64,    // 边长，直径
    blk_cells: usize, // 块单元格数量
}

impl Sender {
    // 构造函数，初始化 Sender 实例
    pub fn new(
        dim: usize,
        radius: u64,
        num_item: u64,
        pk_rec: RistrettoPoint,
        apart: bool,
        sigma: bool,
    ) -> Self {
        let blk_cells: usize = 1 << dim;
        let side_len = 2 * radius;
        let r_l2 = radius * radius;

        let window: usize = if apart == true { 1 } else { blk_cells }; // 根据 apart 参数设置窗口大小
        let m = num_item * window as u64; // 计算消息数量
        let mut rng = rand::thread_rng(); // 创建随机数生成器
        let pk = pk_rec; // 设置公钥
        let mut _coins: Vec<(RistrettoPoint, RistrettoPoint, Scalar)> =
            Vec::with_capacity(m as usize); // 初始化coin向量

        for _ in 0..m {
            // 生成coin
            let a = Scalar::random(&mut rng); // 随机生成一个 Scalar
            let b = Scalar::random(&mut rng); // 随机生成另一个 Scalar

            // _coins.0=a,_coins.1=a*pk,_coins.2=b
            _coins.push((&a * RISTRETTO_BASEPOINT_TABLE, a * pk_rec, b)); // 将coin添加到向量中
        }

        // println!("sender init end");

        return Sender {
            // 返回 Sender 实例
            dim,
            radius,
            sigma,
            m,
            window,
            pk,
            _coins,
            r_l2,
            side_len,
            blk_cells,
        };
    }

    // 发送单个消息
    #[inline]
    pub fn inf_send_msg_single_apart(
        &self,
        encodings: &Vec<okvs::Encoding>, // 编码向量
        pt: &Point,                      // 点
        index: usize,                    // 索引
    ) -> okvs::PointPair {
        let coins = &self._coins[index]; // 获取当前coin
        let mut uv: okvs::PointPair = (RistrettoPoint::identity(), RistrettoPoint::identity()); // 初始化结果
        let mut tem: okvs::PointPair; // 临时变量
        let cel = cell(pt, self.dim, self.radius, self.sigma); // 计算块的左下角坐标
        let key = hash64(&cel); // 计算块的哈希值
                                // 遍历每个维度

        for j in 0..self.dim {
            let key_ij: u64 = hash64(&(pt[j] as u64)); // 计算每个点的哈希值
            tem = okvs::okvs_decode(&encodings[j], key ^ key_ij); // 解码
            uv.0 += tem.0; // 更新结果
            uv.1 += tem.1; // 更新结果
        }
        // 完成计算
        uv.0 = coins.2 * uv.0 + coins.0; // 最终计算
        uv.1 = coins.2 * uv.1 + coins.1; // 最终计算

        return uv; // 返回结果
    }

    // 刷新 Sender 实例
    pub fn refresh(&mut self) {
        let mut rng = rand::thread_rng(); // 创建随机数生成器
        self._coins.clear(); // 清空coin向量
        for _ in 0..self.m {
            // 重新生成coin
            let a = Scalar::random(&mut rng); // 随机生成一个 Scalar
            self._coins.push((
                &a * RISTRETTO_BASEPOINT_TABLE, // 计算coin
                a * self.pk,                    // 计算coin的公钥
                Scalar::random(&mut rng),       // 随机生成一个 Scalar
            ));
        }
    }

    pub fn print(&self) {
        println!(
            "Sender: dim={}, radius={}, sigma={}, m={}, window={}",
            self.dim, self.radius, self.sigma, self.m, self.window
        );
        println!(
            "r_l2: {}, side_len: {}, blk_cells: {}",
            self.r_l2, self.side_len, self.blk_cells
        );
    }
}
