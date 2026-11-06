# S3 dt/mesh 对比隔离 workspace

日期：2026-11-02。用途：在当前生产 `simulator/cpp` 的隔离副本中加入 `dynamic_bar_check`，运行 M2/M3/M4 × dt 的时间离散对比。生产目录未被修改。

## 构建

```bash
cmake -S cpp -B /tmp/fiction-s3-dtmesh-build \
  -DCMAKE_BUILD_TYPE=Release -DFICTION_USE_VIEWER=OFF
cmake --build /tmp/fiction-s3-dtmesh-build --target dynamic_bar_check -j
```

## 运行

```bash
/tmp/fiction-s3-dtmesh-build/dynamic_bar_check \
  --mesh <run>/meshes/M2.msh --out-dir <case-dir> \
  --dt 0.005 --steps 100 --tol 1e-6 --max-iter 50 \
  --mode bend --gravity-z -9.81
```

M3/M4 使用相同命令，仅替换网格、dt 和步数（dt=0.0025 为 200 步，dt=0.00125 为 400 步）。

## 分析

```bash
/home/steven/miniconda3/envs/elastogen/bin/python tools/analyze_s3_dtmesh.py \
  --run simulator/output/validation/s3-dtmesh-20261102-01 \
  --figures figures/validation/s2-s3-refinement-20261102-01
```

## 比较设置

- 在固定空间模型上比较 0–0.5 s 的时间步敏感性。
- M2/M3/M4 分别记录峰值与周期；每个网格内部比较时间步加密。
- 原始 S3 的 2% 逐点对照结果保留为未通过。
