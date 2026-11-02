# 回滚与复现

原工程与第一轮 Variant-A 从未修改，无需对它们回滚。当前 instrumented 副本只增加诊断与参数入口，原算法仍为 Variant-A。

在独立 `rollback-check/instrumented` 副本中已实际逆向应用三文件补丁，并用 cmp 验证三文件均与第一轮 Variant-A 一致。保留 instrumented 供审查。若以后需要撤销该副本的诊断：

```bash
cd /home/steven/mjx/fiction/experiments/workspaces/fiction-s1-two-group-20261031-01/instrumented
patch -R --dry-run -p1 < ../diagnostics.patch
patch -R -p1 < ../diagnostics.patch
```

此逆补丁只恢复到第一轮修正版，不会撤销第一轮“增刚度后重算方向”的修正。请在上述副本目录执行，不要在原工程应用。

复现先使用 CMake 配置和构建 instrumented/cpp 至 `/tmp/fiction-s1-two-group-build`（Release、CHOLMOD），再执行 `run-groups.sh`。新输出写入独立 reproductions 子目录，不覆盖本轮 results。`executed-run-groups.sh` 保留本次执行命令作为记录，不应直接重跑以免覆盖冻结输出。
