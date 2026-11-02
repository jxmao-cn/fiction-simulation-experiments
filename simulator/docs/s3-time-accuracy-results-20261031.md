# S3 短窗口时间精度候选：执行结果

日期：2026-10-31。状态：**reference_not_established**。原S3未通过状态不变。

执行设计见 [候选规格](s3-time-accuracy-revision-20261031.md)。采用原step判据；只更改CSV协议标签，未采用新combined算法。生产源码与旧run未改。

批次数值墙钟 166.82 s，实际完成 2 个病例。每病例180s、总900s上限。

| 病例 | dt | tol | 完成步数 | 退出码 | 数值通过 |
|---|---:|---:|---:|---:|---|
| ref | 0.0003125 | 1e-06 | 1600 | 0 | True |
| refcheck | 0.00015625 | 1e-06 | 3200 | 0 | True |

## 资格与停止

```json
{
  "status": "reference_not_established",
  "results": [
    {
      "case": "ref",
      "dt": 0.0003125,
      "tol": 1e-06,
      "exitcode": 0,
      "steps": 1600,
      "numerical_pass": true
    },
    {
      "case": "refcheck",
      "dt": 0.00015625,
      "tol": 1e-06,
      "exitcode": 0,
      "steps": 3200,
      "numerical_pass": true,
      "reference_check": {
        "max_m": 0.0005750772541847043,
        "rms_m": 0.00023728436077306034
      }
    }
  ],
  "wall_seconds": 166.8193147800048,
  "original_S3_status": "unchanged: failed",
  "production_modified": false
}
```

未执行病例不会记为通过。若数值失败，先停止并保留现场；不放宽tol、不追加迭代追绿。参考或代数资格失败时，不进行粗档精度判定。

证据：[s3-time-accuracy-20261031-01](../output/validation/s3-time-accuracy-20261031-01/)，含spec/config、源码快照、二进制及输入哈希、构建/CTest、实际命令、逐步CSV、末态、退出码及time记录。

本结果仅适用于候选0–0.5s窗口；人工阻尼、相位和振幅差都属于时间误差，未做扣除或相位对齐。
