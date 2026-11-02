#!/usr/bin/env python3
"""生成保形 Freudenthal 六面体分割方杆的 MshIO 4.1 四面体网格。"""

import argparse
import itertools
import json
import math
from collections import Counter
from pathlib import Path


LENGTH = 1.0
WIDTH = 0.05


def determinant(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def make_mesh(nx, ny, nz):
    """在全局一致的体对角线上使用六个轴向排列四面体剖分每个单元。"""
    coords = []
    node_id = {}
    for k in range(nz + 1):
        for j in range(ny + 1):
            for i in range(nx + 1):
                node_id[i, j, k] = len(coords)
                coords.append((LENGTH * i / nx, WIDTH * (j / ny - 0.5), WIDTH * (k / nz - 0.5)))

    permutations = list(itertools.permutations(range(3)))
    tets = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                base = (i, j, k)
                for perm in permutations:
                    p1 = list(base)
                    p1[perm[0]] += 1
                    p2 = p1.copy()
                    p2[perm[1]] += 1
                    p3 = list(base)
                    for d in range(3):
                        p3[d] += 1
                    tet = [node_id[base], node_id[tuple(p1)], node_id[tuple(p2)], node_id[tuple(p3)]]
                    a, b, c = (tuple(coords[tet[m]][d] - coords[tet[0]][d] for d in range(3))
                               for m in (1, 2, 3))
                    if determinant(a, b, c) < 0.0:
                        tet[1], tet[2] = tet[2], tet[1]
                    tets.append(tuple(tet))

    face_orientations = {}
    tet_volumes = []
    for tet in tets:
        p = [coords[v] for v in tet]
        a = tuple(p[1][d] - p[0][d] for d in range(3))
        b = tuple(p[2][d] - p[0][d] for d in range(3))
        c = tuple(p[3][d] - p[0][d] for d in range(3))
        signed = determinant(a, b, c) / 6.0
        tet_volumes.append(signed)
        # 正定向四面体的有向边界，供后续表面渲染保持外法向一致。
        oriented_faces = ((tet[1], tet[2], tet[3]), (tet[0], tet[3], tet[2]),
                          (tet[0], tet[1], tet[3]), (tet[0], tet[2], tet[1]))
        for face in oriented_faces:
            key = tuple(sorted(face))
            perm = [key.index(v) for v in face]
            inversions = sum(perm[i] > perm[j] for i in range(3) for j in range(i + 1, 3))
            sign = -1 if inversions % 2 else 1
            count, first, sign_sum = face_orientations.get(key, (0, face, 0))
            face_orientations[key] = (count + 1, first, sign_sum + sign)

    face_counts = Counter({key: entry[0] for key, entry in face_orientations.items()})
    boundary = [first for count, first, _ in face_orientations.values() if count == 1]
    internal_faces_oppositely_oriented = all(sign_sum == 0 for count, _, sign_sum in face_orientations.values()
                                             if count == 2)
    edge_counts = Counter()
    for face in boundary:
        for q in range(3):
            edge_counts[tuple(sorted((face[q], face[(q + 1) % 3])))] += 1

    total_volume = sum(tet_volumes)
    expected_volume = LENGTH * WIDTH * WIDTH
    expected_boundary_faces = 4 * (nx * ny + nx * nz + ny * nz)
    def on_outer_plane(face):
        pts = [coords[v] for v in face]
        planes = ((0, (0.0, LENGTH)), (1, (-WIDTH / 2, WIDTH / 2)),
                  (2, (-WIDTH / 2, WIDTH / 2)))
        return any(all(abs(p[d] - plane) < 1e-14 for p in pts)
                   for d, values in planes for plane in values)
    checks = {
        "all_positive_signed_volumes": all(v > 0.0 for v in tet_volumes),
        "no_face_shared_by_more_than_two_tets": all(v in (1, 2) for v in face_counts.values()),
        "shared_faces_have_opposite_orientation": internal_faces_oppositely_oriented,
        "closed_boundary_edges": bool(edge_counts) and all(v == 2 for v in edge_counts.values()),
        "boundary_faces_on_outer_box": all(on_outer_plane(face) for face in boundary),
        "boundary_face_count_matches_box": len(boundary) == expected_boundary_faces,
        "volume_matches_box": abs(total_volume - expected_volume) <= 1e-12 * expected_volume,
        "bbox_matches_requested_bar": (
            abs(min(p[0] for p in coords)) < 1e-14
            and abs(max(p[0] for p in coords) - LENGTH) < 1e-14
            and abs(max(p[1] for p in coords) - WIDTH / 2) < 1e-14
            and abs(min(p[1] for p in coords) + WIDTH / 2) < 1e-14
            and abs(max(p[2] for p in coords) - WIDTH / 2) < 1e-14
            and abs(min(p[2] for p in coords) + WIDTH / 2) < 1e-14),
    }
    if not all(checks.values()):
        raise ValueError(f"网格自检失败: {checks}")
    return coords, tets, boundary, {
        "nx": nx, "ny": ny, "nz": nz,
        "nodes": len(coords), "tets": len(tets), "boundary_faces": len(boundary),
        "expected_boundary_faces": expected_boundary_faces,
        "volume": total_volume, "expected_volume": expected_volume,
        "min_signed_tet_volume": min(tet_volumes), "max_face_incidence": max(face_counts.values()),
        "boundary_edge_incidence_values": sorted(set(edge_counts.values())),
        "checks": checks,
    }


def write_msh(path, coords, tets, boundary):
    with path.open("w", encoding="ascii") as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
        f.write(f"$Nodes\n1 {len(coords)} 1 {len(coords)}\n3 0 0 {len(coords)}\n")
        for i in range(len(coords)):
            f.write(f"{i + 1}\n")
        for p in coords:
            f.write(f"{p[0]:.17g} {p[1]:.17g} {p[2]:.17g}\n")
        f.write("$EndNodes\n")
        f.write(f"$Elements\n1 {len(tets)} 1 {len(tets)}\n3 0 4 {len(tets)}\n")
        for i, tet in enumerate(tets):
            f.write(f"{i + 1} {tet[0] + 1} {tet[1] + 1} {tet[2] + 1} {tet[3] + 1}\n")
        f.write("$EndElements\n")
        f.write(f"$Surface\n{len(boundary)}\n")
        for tri in boundary:
            f.write(f"{tri[0] + 1} {tri[1] + 1} {tri[2] + 1}\n")
        f.write("$EndSurface\n")


def main():
    parser = argparse.ArgumentParser(description="生成 S1 结构化方杆 MshIO 4.1 tet 网格")
    parser.add_argument("--nx", type=int, required=True, help="x 方向单元数")
    parser.add_argument("--ny", type=int, required=True, help="y 方向单元数")
    parser.add_argument("--nz", type=int, required=True, help="z 方向单元数")
    parser.add_argument("--output", required=True, help="输出 .msh 路径")
    parser.add_argument("--metadata", help="可选自检 JSON 输出路径")
    args = parser.parse_args()
    if min(args.nx, args.ny, args.nz) <= 0:
        parser.error("nx、ny、nz 必须为正整数")
    coords, tets, boundary, report = make_mesh(args.nx, args.ny, args.nz)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_msh(output, coords, tets, boundary)
    if args.metadata:
        metadata = Path(args.metadata)
        metadata.parent.mkdir(parents=True, exist_ok=True)
        metadata.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"[OK] {output}: nodes={len(coords)} tets={len(tets)} boundary_faces={len(boundary)} "
          f"volume={report['volume']:.12g}; checks=PASS")


if __name__ == "__main__":
    main()
