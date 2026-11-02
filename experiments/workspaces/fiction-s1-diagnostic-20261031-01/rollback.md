# 回滚说明

原工程 `/home/steven/mjx/fiction/simulator` 从未修改，因此不需要对原工程回滚，也不要把补丁合并回原工程。

Variant-A 保留修改供审查。若审查后要在 Variant-A 撤销 A，先执行逆向 dry-run，再实际逆向应用：

```bash
cd /home/steven/mjx/fiction/experiments/workspaces/fiction-s1-diagnostic-20261031-01/variant-a
patch -R --dry-run -p1 < ../fix-a.patch
patch -R -p1 < ../fix-a.patch
sha256sum cpp/src/time_integrator.cpp
```

撤销后预期 SHA-256：`84ee9de8ed6118c748d5625b1ddbea7f4b632cb909ac7f8c22978553a9216dc6`，与 `baseline/cpp/src/time_integrator.cpp` 相同。本轮已在独立 `rollback-check/` 副本实际执行这组逆向验证；Variant-A 本身没有撤销，以保留可审查状态。不要删除整个诊断目录。
