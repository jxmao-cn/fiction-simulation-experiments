#!/usr/bin/env python3
# =============================================================================
# tools/bgeo_to_msh.py — Elastogen 体素网格 (bgeo) → MshIO 4.1 ASCII 四面体网格
# =============================================================================
# 用途：
#   把 Elastogen 论文的 cantilever 杆网格（neuralProj/assets/rod*.bgeo，六面体体素网格）
#   转成 fiction 可直接加载的 .msh（MshIO 4.1 ASCII，tet4），用于复现
#   {rods.pdf 重力弯曲} 与 {rods_twist.pdf 扭转} 实验（见 docs/elastogen-rods-experiments.md）。
#
# 转换规则：
#   - 读取 bgeo 的六面体拓扑（pyHouGeoIO.GetCubeTopology -> (8, nHex)）
#   - 按局部 bit 顺序识别每个六面体 8 角（bit0=x, bit1=y, bit2=z）
#   - Kuhn(Freudenthal) 6-tet 剖分：全体共享体对角线，按 (i+j+k) 奇偶翻转主对角，
#     保证相邻六面体共享面上三角对齐（保形网格）
#   - 顶点按坐标字典序稳定重编号，tag 1-based
#
# 自检（默认失败即不写盘）：
#   1. 每个体素为 2x2x2 轴对齐角格；8 顶点齐全
#   2. 每六面体的 6-tet 体积和 == 六面体体积（相对容差 1e-9）
#   3. 无退化四面体（|V_tet| > 1e-12）
#   4. 边界三角面闭合：每条无向边恰好被 2 个边界面共享
#
# 运行（需 elastogen conda 环境以匹配 pyHouGeoIO.cp310 扩展）：
#   conda run -n elastogen python tools/bgeo_to_msh.py \
#       --input  ../Elastogen/neuralProj/assets/rod.bgeo \
#       --output cpp/tests/data/rod16x3x3.msh
# =============================================================================

import argparse
import os
import sys
from collections import defaultdict

import numpy as np

# Kuhn 6-tet（偶奇偶基准）：全部共享主对角线 0↔7
# 局部位序: bit0=x, bit1=y, bit2=z → 0:(0,0,0) 1:(1,0,0) 2:(0,1,0) 3:(1,1,0)
#                    4:(0,0,1) 5:(1,0,1) 6:(0,1,1) 7:(1,1,1)
KUHN_EVEN = (
    (0, 7, 1, 3),
    (0, 7, 3, 2),
    (0, 7, 2, 6),
    (0, 7, 6, 4),
    (0, 7, 4, 5),
    (0, 7, 5, 1),
)
# 奇偶性翻转：bit 全取反（镜像六面体），使跨体共享面剖分一致
BIT_FLIP = tuple(7 - b for b in range(8))


def load_bgeo(path, pyhgeoio_dir):
    if pyhgeoio_dir and pyhgeoio_dir not in sys.path:
        sys.path.insert(0, pyhgeoio_dir)
    import pyHouGeoIO

    cube = pyHouGeoIO.HCube()
    pyHouGeoIO.ImportHouGeo(path, cube)
    P = np.asarray(cube.GetPointAttributeMatXf('P'))
    topo = np.asarray(cube.GetCubeTopology())
    if topo.shape[0] != 8:
        topo = topo.T
    return P, topo


def tet_volume(P, tet):
    a, b, c, d = (P[i] for i in tet)
    return abs(float(np.dot(b - a, np.cross(c - a, d - a)))) / 6.0


def build_boundary(tets):
    """边界三角面（四面体面只出现一次的），返回定向三角面列表（排序稳定）"""
    count = defaultdict(int)
    oriented = {}
    faces = ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))
    for tet in tets:
        for f in faces:
            tri = tuple(tet[k] for k in f)
            key = tuple(sorted(tri))
            count[key] += 1
            oriented[key] = tri
    return sorted(oriented[k] for k, n in count.items() if n == 1)


def convert(P, topo):
    """返回 (Ps 排序后坐标, tets 0-based, boundary, errors)"""
    errors = []
    n = P.shape[0]
    order = np.lexsort((P[:, 2], P[:, 1], P[:, 0]))
    remap = np.empty(n, dtype=np.int64)
    remap[order] = np.arange(n)
    Ps = P[order]

    tets = []
    for hi in range(topo.shape[1]):
        verts = topo[:, hi]
        pts = P[verts]
        xs, ys, zs = (np.unique(pts[:, k]) for k in range(3))
        if not (len(xs) == len(ys) == len(zs) == 2):
            errors.append(f"hex {hi}: 非轴对齐 2x2x2 角格")
            continue
        mn = np.array([xs[0], ys[0], zs[0]])
        mx = np.array([xs[1], ys[1], zs[1]])
        bits = {}
        for v in verts:
            d = (P[v] - mn) / (mx - mn)
            key = tuple(int(np.rint(d[b])) for b in range(3))
            bits[key] = int(v)
        if len(bits) != 8:
            errors.append(f"hex {hi}: 角格不满 8 顶点")
            continue
        corner_ids = [bits[(b & 1, (b >> 1) & 1, (b >> 2) & 1)] for b in range(8)]

        ijk = tuple(int(round(c)) for c in mn)
        parity = (ijk[0] + ijk[1] + ijk[2]) % 2
        cell_vol = float(np.prod(mx - mn))
        sum_vol = 0.0
        for a, b, c, d in KUHN_EVEN:
            if parity:
                a, b, c, d = BIT_FLIP[a], BIT_FLIP[b], BIT_FLIP[c], BIT_FLIP[d]
            tet = tuple(int(remap[corner_ids[k]]) for k in (a, b, c, d))
            v = tet_volume(Ps, tet)
            if v < 1e-12:
                errors.append(f"hex {hi}: 退化四面体 vol={v}")
            sum_vol += v
            tets.append(tet)
        if abs(sum_vol - cell_vol) > 1e-9 * max(1.0, cell_vol):
            errors.append(f"hex {hi}: 6-tet 体积和 {sum_vol:.12g} != 六面体 {cell_vol:.12g}")

    boundary = build_boundary(tets)
    edges = defaultdict(int)
    for f in boundary:
        for k in range(3):
            a, b = sorted((f[k], f[(k + 1) % 3]))
            edges[(a, b)] += 1
    if not all(c == 2 for c in edges.values()):
        bad = sum(1 for c in edges.values() if c != 2)
        errors.append(f"边界提取不闭合（{bad} 条边被 !=2 个面共享）")

    return Ps, np.asarray(tets, dtype=np.int64), boundary, errors


def write_msh(path, Ps, tets, boundary):
    n = Ps.shape[0]
    with open(path, "w") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n1 {n} 1 {n}\n3 0 0 {n}\n")
        for i in range(n):
            f.write(f"{i + 1}\n")
        for x, y, z in Ps:
            f.write(f"{x:.9g} {y:.9g} {z:.9g}\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n1 {len(tets)} 1 {len(tets)}\n3 0 4 {len(tets)}\n")
        for e, t in enumerate(tets):
            f.write(f"{e + 1} {t[0] + 1} {t[1] + 1} {t[2] + 1} {t[3] + 1}\n")
        f.write("$EndElements\n")
        f.write(f"$Surface\n{len(boundary)}\n")
        for tri in boundary:
            f.write(f"{tri[0] + 1} {tri[1] + 1} {tri[2] + 1}\n")
        f.write("$EndSurface\n")


def main():
    ap = argparse.ArgumentParser(description="bgeo(HCube) -> MshIO 4.1 ASCII tet 网格")
    ap.add_argument("--input", required=True, help="输入 .bgeo（六面体体素网格）")
    ap.add_argument("--output", required=True, help="输出 .msh 路径")
    ap.add_argument("--pyhgeoio-dir", default="/home/steven/mjx/Archive/Elastogen/neuralProj/utils",
                    help="pyHouGeoIO 扩展模块所在目录")
    ap.add_argument("--allow-errors", action="store_true", help="存在校验告警时仍写盘（默认不写）")
    args = ap.parse_args()

    P, topo = load_bgeo(args.input, args.pyhgeoio_dir)
    print(f"[bgeo] {args.input}: {P.shape[0]} 顶点, {topo.shape[1]} 六面体")

    Ps, tets, boundary, errors = convert(P, topo)
    if errors:
        for e in errors[:10]:
            print("  [ERROR]", e)
        if not args.allow_errors:
            raise SystemExit("转换自检失败，未写盘（--allow-errors 可强制写出）")

    total_vol = sum(tet_volume(Ps, tuple(int(x) for x in t)) for t in tets)
    bbox_min, bbox_max = Ps.min(axis=0), Ps.max(axis=0)
    bbox_vol = float(np.prod(bbox_max - bbox_min))
    print(f"  [OK] tets={len(tets)}, 边界面={len(boundary)}")
    print(f"  [OK] 总体积={total_vol:.12g}, bbox体积={bbox_vol:.12g}, bbox={tuple(bbox_min)}-{tuple(bbox_max)}")
    if abs(total_vol - bbox_vol) > 1e-9 * bbox_vol:
        raise SystemExit(f"总体积 {total_vol} != 包围盒 {bbox_vol}，请检查")

    write_msh(args.output, Ps, tets, boundary)
    print(f"  [OK] 写出 {args.output}（{Ps.shape[0]} 顶点 / {len(tets)} 四面体 / {len(boundary)} 边界面）")


if __name__ == "__main__":
    main()
