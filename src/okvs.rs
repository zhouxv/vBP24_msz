const EPSILON: f64 = 0.5;
const STAT_BITS: u64 = 40;
const FACTOR: f64 = STAT_BITS as f64 * 1.44 as f64;
const HASH_SEED: u64 = 0x1234567890abcdef;
const KAPPA: u64 = 40;

use curve25519_dalek::constants::RISTRETTO_BASEPOINT_TABLE;
use curve25519_dalek::traits::Identity;
use curve25519_dalek::RistrettoPoint;
use curve25519_dalek::Scalar;

// use counting_sort::CountingSort;
use fxhash::hash64;
use rand::rngs::OsRng;

// 键值对
pub type PointPair = (RistrettoPoint, RistrettoPoint);
// OKVS 输入
pub type Encoding = Vec<PointPair>;

// 排序算法
pub fn counting_sort<F, T>(arr: &mut Vec<T>, min: usize, max: usize, key: F)
where
    F: Fn(&T) -> usize,
    T: Clone,
{
    let mut prefix_sums = {
        // 1. Initialize the count array with default value 0.
        let len = max - min;
        let mut count_arr = Vec::with_capacity(len);
        count_arr.resize(len, 0);

        // 2. Scan elements to collect counts.
        for value in arr.iter() {
            count_arr[key(value)] += 1;
        }

        // 3. Calculate prefix sum.
        count_arr
            .into_iter()
            .scan(0, |state, x| {
                *state += x;
                Some(*state - x)
            })
            .collect::<Vec<usize>>()
    };

    // 4. Use prefix sum as index position of output element.
    for value in arr.to_vec().iter() {
        let index = key(value);
        arr[prefix_sums[index]] = value.clone();
        prefix_sums[index] += 1;
    }
}

pub struct GBF {
    data: Vec<RistrettoPoint>,
    m: u64,
    n: u64,
}

impl GBF {
    pub fn new(num_item: u64) -> Self {
        let m: u64 = (num_item as f64 * FACTOR).floor() as u64;
        let data: Vec<RistrettoPoint> =
            (0..m).map(|_| RistrettoPoint::random(&mut OsRng)).collect();
        return GBF {
            data,
            m,
            n: num_item,
        };
    }
    pub fn num_items(&self) -> u64 {
        return self.n;
    }
    pub fn len(&self) -> u64 {
        return self.m;
    }
    pub fn num_hashes(&self) -> u64 {
        return KAPPA;
    }
    pub fn encode(&mut self, list: &Vec<(u64, RistrettoPoint)>) {
        let mut touched: Vec<bool> = vec![false; self.m as usize];
        let mut hashkey: [u64; STAT_BITS as usize] = [0u64; STAT_BITS as usize];
        for (key, value) in list.iter() {
            let seed = hash64(&(*key ^ HASH_SEED));
            for i in 0..STAT_BITS {
                hashkey[i as usize] = hash64(&(seed ^ i)) % self.m;
            }
            let mut temp = *value;
            let mut pos = 0usize;
            for i in 0..STAT_BITS as usize {
                let j = hashkey[i] as usize;
                if false == touched[j] {
                    pos = j;
                    touched[j] = true;
                }
                temp -= self.data[j];
            }
            assert!(pos != 0usize);
            self.data[pos] += temp;
        }
    }

    pub fn decode(&mut self, key: u64) -> RistrettoPoint {
        let seed: u64 = hash64(&(key ^ HASH_SEED));
        let mut hashkey = [0u64; STAT_BITS as usize];
        for i in 0..STAT_BITS {
            hashkey[i as usize] = hash64(&(seed ^ i)) % self.m;
        }

        let mut result = RistrettoPoint::identity();
        for i in 0..STAT_BITS as usize {
            result += self.data[hashkey[i] as usize];
        }
        return result;
    }
}

pub fn okvs_decode(data: &Encoding, key: u64) -> PointPair {
    let pos_band_range: u64 = data.len() as u64 - KAPPA;
    let seed: u64 = hash64(&(key ^ HASH_SEED));
    let pos: usize = (seed % pos_band_range) as usize;
    let band: Vec<u64> = (0..KAPPA).map(|i| (seed >> i) & 0x01).collect();

    let mut result = (RistrettoPoint::identity(), RistrettoPoint::identity());

    for i in pos..KAPPA as usize + pos {
        if 0 == band[i - pos] {
            continue;
        }
        result.0 += data[i].0;
        result.1 += data[i].1;
    }
    return result;
}

pub fn okvs_decode_batch(data: &Encoding, keys: &Vec<u64>) -> Encoding {
    let pos_band_range: u64 = data.len() as u64 - KAPPA;
    let mut result: Encoding = Vec::with_capacity(keys.len());
    for key in keys {
        // 身份元素,或可叫零点，无穷远点
        let mut result_p = (RistrettoPoint::identity(), RistrettoPoint::identity());
        let seed: u64 = hash64(&(key ^ HASH_SEED));
        let pos: usize = (seed % pos_band_range) as usize;
        let band: Vec<u64> = (0..KAPPA).map(|i| (seed >> i) & 0x01).collect();
        for i in pos..KAPPA as usize + pos {
            if 0 == band[i - pos] {
                continue;
            }
            result_p.0 += data[i].0;
            result_p.1 += data[i].1;
        }
        result.push(result_p);
    }
    return result;
}
pub struct OkvsGen {
    _data: Vec<(Scalar, Scalar)>,
    // Eoncoing 用的矩阵，参数分别为 哈希带起始位置，主元的位置，哈希带，Value
    _matrix: Vec<(usize, usize, (Vec<Scalar>, (Scalar, Scalar)))>,
    // 存储总的 m 值，它是编码的目标大小，通过输入项目数 num_item 计算得出
    m: u64,
    // 存储输入项的数量
    n: u64,
    // 存储输入项的数量
    epsilon: f64,
}
//
impl OkvsGen {
    pub fn new(num_item: u64) -> Self {
        // 根据 n 和 epsilon 计算 m
        let mut m: u64 = (num_item as f64 * (1.0 + EPSILON)).ceil() as u64;

        if m < KAPPA {
            m = (KAPPA as f64 * (1.0 + EPSILON)).ceil() as u64;
        }

        // 随机生成 m 个随机 (Scalar, Scalar)
        let _data: Vec<(Scalar, Scalar)> = (0..m)
            // .map(|_| (Scalar::ZERO, Scalar::ZERO))
            .map(|_| (Scalar::random(&mut OsRng), Scalar::random(&mut OsRng)))
            .collect();
        return OkvsGen {
            _data,
            // _matrix 容量初始化为 n
            _matrix: Vec::with_capacity(num_item as usize),
            m,
            n: num_item,
            epsilon: EPSILON,
        };
    }
    pub fn num_items(&self) -> u64 {
        return self.n;
    }
    pub fn len(&self) -> u64 {
        return self.m;
    }
    pub fn expansion_rate(&self) -> f64 {
        return self.epsilon;
    }

    /*
    清空 _matrix 和 _data，并重新填充 _data。这使得 OkvsGen 可以在每次调用时重新生成新的随机数据
    */
    pub fn refresh(&mut self) {
        self._matrix.clear();
        self._data.clear();
        for _ in 0..self.m {
            self._data
                .push((Scalar::random(&mut OsRng), Scalar::random(&mut OsRng)));
        }
    }

    pub fn encode(&mut self, list: &Vec<(u64, (Scalar, Scalar))>) -> Encoding {
        /*
        初始化 encoding，size为 m, 每个元素 Value 为 Ristretto group for Curve25519 中的两个点
        */
        let mut data: Encoding = Vec::with_capacity(self.m as usize);

        /*
        基于输入的键值对，通过哈希操作生成位置和哈希带，并将结果添加到一个矩阵中
        */
        let pos_band_range: u64 = self.m - KAPPA;
        for (key, value) in list.iter() {
            let seed = hash64(&(*key ^ HASH_SEED));
            let hash_pos = (seed % pos_band_range) as usize;
            // 哈希带长度为 kappa ，每个是一个 Scalar
            let hash_band: Vec<Scalar> = (0..KAPPA)
                .map(|i: u64| Scalar::from((seed >> i) & 0x01))
                .collect();

            self._matrix.push((hash_pos, 0usize, (hash_band, *value)));
        }

        /*
        _matrix 按照哈希位置进行排序，以确保数据按位置有序
        */
        self._matrix.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        // counting_sort(&mut self._matrix, 0usize, pos_band_range as usize, |t| t.0);

        let mut pivots: usize = 0;

        /*
        高斯消元法（主元化和行归一化）
        */
        for row in 0..self.n as usize {
            // 将 _matrix 分割为两部分：top 包含当前行及其之前的所有行，bot 包含当前行之后的所有行
            let (top, bot) = self._matrix.split_at_mut(row + 1);

            /*
            获取top矩阵的最后一行，寻找主元并进行归一化，并对后续行进行处理
             */
            if let Some((pos, piv, (band, value))) = top.last_mut() {
                // 对哈希带的每一位进行处理
                for i in 0..KAPPA as usize {
                    // 如果当前 band 的元素为零，则跳过该元素
                    if Scalar::ZERO == band[i] {
                        continue;
                    }
                    // 更新主元的位置，并增加主元计数
                    *piv = *pos + i;
                    pivots += 1;

                    // if the pivot is not one, divide the row to normalize it
                    // 归一化主元, 如果主元不为1, 要对整行进行归一化
                    if Scalar::ONE != band[i] {
                        let inv = band[i].invert();
                        band[i] = Scalar::ONE;
                        for j in i + 1..KAPPA as usize {
                            band[j] *= inv;
                        }
                        value.0 *= inv;
                        value.1 *= inv;
                    }

                    // update the band from the following rows
                    // for j in row + 1..self.n as usize {
                    for (pos_j, piv_j, (band_j, value_j)) in bot.iter_mut() {
                        // skip if the position is out of range
                        if *pos_j > *piv {
                            break;
                        }
                        // if found the non-zero position, subtract the band
                        *piv_j = *piv - *pos_j;
                        let multplier = band_j[*piv_j];
                        if Scalar::ZERO != multplier {
                            if Scalar::ONE != multplier {
                                for k in i + 1..KAPPA as usize {
                                    band_j[*piv_j + k - i] -= band[k] * multplier;
                                }
                                value_j.0 -= value.0 * multplier;
                                value_j.1 -= value.1 * multplier;
                            } else {
                                for k in i + 1..KAPPA as usize {
                                    band_j[*piv_j + k - i] -= band[k];
                                }
                                value_j.0 -= value.0;
                                value_j.1 -= value.1;
                            }
                        }
                        band_j[*piv_j] = Scalar::ZERO;
                    }
                    break;
                }
            }
        }

        // band should not be all-zero
        // pivots不应该全为0,否则panic
        assert_eq!(pivots, self.n as usize);

        // back substitution，回代
        for (pos, piv, (band, value)) in self._matrix.iter().rev() {
            let (mut val, mut val2) = *value;
            for i in 0..KAPPA as usize {
                if Scalar::ZERO == band[i] {
                    continue;
                }
                // 主元列不需要做更新，因为主元列已经用于消去其它列的元素
                if (i + *pos) == *piv {
                    continue;
                }
                // 主元对应的索引 piv 之后的 已经随机化了，计算主元对应的Encoding的value
                val -= self._data[i + *pos].0 * band[i];
                val2 -= self._data[i + *pos].1 * band[i];
            }
            self._data[*piv] = (val, val2);
        }

        /*
        finalize
        将Scalar转换为RistrettoPoint
        */
        for i in 0..self.m as usize {
            data.push((
                &self._data[i].0 * RISTRETTO_BASEPOINT_TABLE,
                &self._data[i].1 * RISTRETTO_BASEPOINT_TABLE,
            ));
        }

        // 清理 matrix
        self._matrix.clear();
        // self._data.clear();
        return data;
    }
}

#[cfg(test)]
mod tests {
    use crate::okvs;
    use curve25519_dalek::scalar::Scalar;
    use fxhash::hash64;
    use rand::rngs::OsRng; // Add this import for bincode
    use std::time::Instant;

    #[test]
    fn bench_okvs_decode() {
        let n = 83886080;
        let mut list: Vec<(u64, (Scalar, Scalar))> = Vec::new();
        eprintln!("{} items, OkvsGen.Decode", n);
        for j in 0..n {
            list.push((
                hash64(&j),
                (Scalar::random(&mut OsRng), Scalar::random(&mut OsRng)),
            ));
        }

        eprintln!("OKVS开始!!");
        let start1 = Instant::now();
        let mut okvs_instance: okvs::OkvsGen = okvs::OkvsGen::new(n);
        let duration1 = start1.elapsed();
        let data = okvs_instance.encode(&list);
        eprintln!("time1:{}", duration1.as_secs());

        eprintln!("OKVS encoding 结束!!");

        let start2 = Instant::now();
        let keys: Vec<u64> = (0..n).collect();
        okvs::okvs_decode_batch(&data, &keys);
        let duration2 = start2.elapsed();
        eprintln!("time2:{}", duration2.as_secs());

        eprintln!("OKVS decoding 结束!!");
    }
}
