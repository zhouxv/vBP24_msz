## build img

sudo docker build -t vbp24_upsi:latest .
docker tag vbp24_upsi:latest blueobsidian/vbp24_upsi:latest
docker push blueobsidian/vbp24_upsi:latest

## run container

sudo docker run -dit --name vbp24_100Mbps --cap-add=NET_ADMIN --cpuset-cpus="0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30" --cpus=16 --memory=128g --memory-swap=128g vbp24_upsi:latest
sudo docker run -dit --name vbp24_10Mbps --cap-add=NET_ADMIN --cpuset-cpus="32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62" --cpus=16 --memory=128g --memory-swap=128g vbp24_upsi:latest

## tcconfig

In our image, we use [tcconfig](!https://github.com/thombashi/tcconfig) to set up traffic control of network bandwidth/latency. Usage instructions are provided below.

```bash
tcset lo --rate 10Gbps --overwrite               # Set the local loopback interface bandwidth to 10Gbps
tcset lo --rate 1Gbps --delay 5ms --overwrite    # Set bandwidth to 1Gbps and add 5ms network delay
tcset lo --rate 100Mbps --delay 80ms --overwrite # Set bandwidth to 100Mbps and add 20ms network delay
tcset lo --rate 10Mbps --delay 80ms --overwrite # Set bandwidth to 100Mbps and add 20ms network delay
tcshow lo                                        # Display current traffic control settings for the loopback interface
tcdel lo -a                                      # Remove all traffic control rules from the loopback interface
```

```
nohup ./run_bench.sh > vbp24_100Mbps.log 2>&1 &
```
