FROM ubuntu:22.04

RUN apt-get update
RUN apt-get install -y net-tools iproute2 python3 python3-pip
RUN pip install tcconfig
RUN apt-get install -y build-essential curl

WORKDIR /home

RUN export RUSTUP_DIST_SERVER=https://mirrors.ustc.edu.cn/rust-static && \
    export RUSTUP_UPDATE_ROOT=https://mirrors.ustc.edu.cn/rust-static/rustup && \
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

COPY ./src /home/src
COPY ./Cargo.toml /home
COPY ./*.sh /home


RUN chmod +x /home/*.sh && \
    export PATH="$HOME/.cargo/bin:$PATH" && \
    cargo build && \
    cp /home/target/debug/vBP24_ufpsi /home/vBP24_ufpsi





