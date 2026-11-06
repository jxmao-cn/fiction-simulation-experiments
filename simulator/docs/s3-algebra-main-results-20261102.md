# S3代数求解精度资格：主工程执行结果

日期2026-11-02。状态：**numerical_or_resource_failure**。原S3未通过与旧参考未建立状态保持。

采用真实生产core的step模式，材料/加载/几何固定；运行前规格见 [规格](s3-algebra-main-spec-20261102.md)。

```json
{
  "status": "numerical_or_resource_failure",
  "results": [
    {
      "case": "ref",
      "dt": 0.0003125,
      "tol": 1e-06,
      "steps": 1600,
      "exitcode": 0,
      "numerical_pass": true
    },
    {
      "case": "ref-tight",
      "dt": 0.0003125,
      "tol": 1e-07,
      "steps": 20,
      "exitcode": 1,
      "numerical_pass": false
    }
  ],
  "wall_seconds": 59.08617666497594,
  "original_S3_status": "unchanged: failed",
  "production_modified": false
}
```

松紧容差门槛：101共同物理时刻尖端最大差与完整末态节点最大差均<=10 μm。数值、资源或门槛失败即停。两时间步差用于诊断；本批检查代数求解精度。

证据：[运行目录](../output/validation/s3-algebra-main-20261102-01/)，保存源码快照、哈希、生产core链接/头依赖、构建/CTest、命令/CSV/末态/资源/逐采样差值。未改生产代码。

## 首次失败现场

紧容差病例在第20步（t=0.00625 s）达到50次Newton上限，方向率3.591627476e-7高于tol=1e-7，累计线搜索减半1546次；lineSearchFailed=false、状态有限、DBC误差0、min J=0.999817854。终止原因为迭代门槛。raw IP梯度1.260525869e-11，换算物理残差约1.29078e-4 N。停滞机制及combined对照留待逐轮trace诊断。

未运行dt=0.00015625的两病例，未改变maxIter、tol、模式或重跑；由于tight只得到0.00625s前缀，无法完成0–0.5s松紧容差资格比较。因此结论为代数求解精度资格未建立，原参考未建立仍保留。进一步诊断须另立批次，本轮停止。

逐case原runner的JSON/CSV可能保留旧protocol/spec标签（观测代码原样复用）；本次执行身份以根目录config.json、spec.md、实际命令、源码快照和本文为准。
