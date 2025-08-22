use std::u128;

use rand::Rng;
use vbp24_ufpsi::protocol;
use vbp24_ufpsi::psi::Point;

use clap::Parser;
#[derive(Parser)]
#[command(version, about, long_about = None)]
struct Cli {
    /// 维度
    #[arg(short, long, default_value_t = 2)]
    dim: usize,

    /// 半径
    #[arg(short, long, default_value_t = 10)]
    radius: usize,

    /// sender方集合大小 log
    #[arg(short, long, default_value_t = 8)]
    n: usize,

    /// receiver方集合大小 log
    #[arg(short, long, default_value_t = 12)]
    m: usize,

    /// 运行次数
    #[arg(short, long, default_value_t = 1)]
    times: usize,

    /// 2σ or 4σ, 默认2σ
    #[arg(short, long)]
    sigma: bool,

    #[arg(short, long, default_value_t = 8080)]
    port: usize,
}

fn main() {
    let cli = Cli::parse();

    let dim = cli.dim;
    let radius = cli.radius;
    let n = cli.n;
    let m = cli.m;
    let times = cli.times;
    let port = cli.port;
    let sigma = cli.sigma;

    println!(
        "dim: {}, radius: {}, n:{} {}, m:{} {}, times: {}, port: {}, sigma: {}",
        dim,
        radius,
        n,
        2usize.pow(n as u32),
        m,
        2usize.pow(m as u32),
        times,
        port,
        sigma
    );

    lin_low(
        dim,
        radius,
        2usize.pow(n as u32),
        2usize.pow(m as u32),
        times,
        port,
        sigma,
    );
}

fn lin_low(
    dim: usize,
    radius: usize,
    n_s: usize,
    n_r: usize,
    times: usize,
    port: usize,
    sigma: bool,
) {
    let mut a = 0;
    let mut b = 0;
    let mut c = 0;

    for _ in 0..times {
        let mut data_s = sample_test_data_points(n_s, dim, radius);
        let data_r = sample_test_data_points(n_r, dim, radius);

        for i in 0..dim {
            data_s[9][i] = data_r[7][i] - (radius / 2) as u64;
            data_s[11][i] = data_r[7][i] + radius as u64;
        }

        let (psi_rcr, psi_sdr) = protocol::setup(dim, radius, n_s, n_r, true, sigma);

        let (d, e, f) = protocol::run_standard(psi_rcr, psi_sdr, data_r, data_s, port);

        // println!("msg1 {} msg2 {} times {}", d, e, f);
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

fn sample_test_data_points(num: usize, dim: usize, radius: usize) -> Vec<Point> {
    let side_len = 2 * radius as u64;
    let mut rng = rand::thread_rng();
    let mut points: Vec<Point> = Vec::with_capacity(num);
    for _ in 0..num {
        let mut point: Point = vec![0u64; dim];
        for i in 0..dim {
            point[i] = rng.gen_range(side_len..=(1 << 31));
        }
        points.push(point);
    }
    return points;
}
