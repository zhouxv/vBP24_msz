#!/bin/bash

# 参数数组
ns=(5 8)
ms=(18 20)
rs=(10 30)
ds=(2 6)


# 循环执行
for d in "${ds[@]}"; do
  for n in "${ns[@]}"; do
    for m in "${ms[@]}"; do
      for r in "${rs[@]}"; do
        ./vBP24_ufpsi -n $n -m $m -d $d -r $r 
        echo   # 输出空行
      done
    done
  done
done

# 循环执行
for d in "${ds[@]}"; do
  for n in "${ns[@]}"; do
    for m in "${ms[@]}"; do
      for r in "${rs[@]}"; do
        ./vBP24_ufpsi -n $n -m $m -d $d -r $r -s
        echo   # 输出空行
      done
    done
  done
done