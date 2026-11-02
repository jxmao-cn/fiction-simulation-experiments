#!/usr/bin/env python3
"""make_contact_sheet.py — 从实验视频抽帧拼接带标注的 contact sheet（工作报告第三部分供图）。

对每个实验系列（ipc_models / elastogen_cantilever）：
  * 用 ffmpeg 从各视频抽取起始/中段/末态三张关键帧；
  * 按「行 = 实验、列 = 时间点」拼图，行标签带实验名 + 关键参数（材料 spring/nh、载荷、步数），
    每格左上角再打参数角标；
  * 帧数不齐时裁齐到最短有效帧并标注，绝不留纯黑格（缺帧画浅灰占位 + 文字说明）。

用法（在 fiction/ 仓库根目录）：
    python3 tools/make_contact_sheet.py --series ipc_models
    python3 tools/make_contact_sheet.py --series elastogen_cantilever
    python3 tools/make_contact_sheet.py --all

输出：
    output/videos/<series>/<series>_sheet.png

依赖：系统 python3 自带 PIL + ffmpeg（禁止系统 pip 安装；若 PIL 不可用则报错退出并提示
改用 conda run -n solid-sim python）。视频内嵌参数以 MANIFEST 为准，与
docs/elastogen-rods-experiments.md:50-78、AGENTS.md:91-135 的复现配方一一对应。
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import date

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

try:
    from PIL import Image, ImageDraw, ImageFont
except ImportError:
    sys.exit("[make_contact_sheet] 缺少 PIL；请勿系统 pip 安装，改用 conda run -n solid-sim python 运行")

# ---------------------------------------------------------------- 画布参数
MARGIN = 28          # 外边距
GUTTER = 10          # 格间距
LABEL_W = 420        # 行标签栏宽度
CELL_W, CELL_H = 640, 360   # 单格（16:9，对应 1920x1080 视频帧等比缩小）
HEADER_H = 108       # 顶部标题 + 列头高度
TITLE_H = 58

BG = (18, 21, 28)
PANEL = (30, 35, 46)
PANEL_EDGE = (60, 70, 88)
TEXT = (232, 236, 242)
TEXT_DIM = (150, 160, 174)
ACCENT = (94, 168, 255)
BADGE_BG = (14, 17, 23)

FONT_CJK_CANDIDATES = [
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
    "/usr/share/fonts/opentype/noto/NotoSansCJK-Bold.ttc",
    "/usr/share/fonts/truetype/arphic/ukai.ttc",
]
FONT_MONO = "/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf"


def load_font(size, bold=False):
    paths = FONT_CJK_CANDIDATES[::-1] if bold else FONT_CJK_CANDIDATES
    for p in paths:
        if os.path.exists(p):
            for idx in range(0, 5):
                try:
                    f = ImageFont.truetype(p, size, index=idx)
                    # 中文可渲染性探测
                    if f.getbbox("材料载荷步数")[2] > 0:
                        return f
                except Exception:
                    continue
    return ImageFont.load_default()


# ---------------------------------------------------------------- 实验清单
# 材料/载荷/步数等参数与 CLI 配方一致（禁改；如场景变化请同步 docs 与 README）。
MANIFEST = {
    "ipc_models": {
        "title": "IPC 几何 + fiction 物理 演示系列（ipc_models）",
        "series_note": "网格来自 IPC/input/tetMeshes（SIGGRAPH 2020 资产），物理为 fiction 质量-弹簧 + IP 投影牛顿；960x720 旧版见 output/_archive_960x720/",
        "rows": [
            {
                "file": "gravity_bar.mp4",
                "name": "gravity_bar  重力下垂",
                "material": "spring（边弹簧，k=youngs·l_avg=7.75e5）",
                "material_tag": "spring",
                "load": "重力 g=(0,0,-9.81)，bar-186 左端 4 节点固定",
                "load_tag": "重力 -9.81",
                "steps": "1000 步 × dt=0.01 s（tol 1e-6）",
                "steps_tag": "1000 步",
            },
            {
                "file": "force_bar.mp4",
                "name": "force_bar  区域外力拉伸",
                "material": "spring（边弹簧，k=youngs·l_avg=7.75e5）",
                "material_tag": "spring",
                "load": "区域外力 a=(0,0,40) 作用于右端盒（IPC NBC 等价），无重力",
                "load_tag": "外力 (0,0,40)",
                "steps": "1000 步 × dt=0.01 s（tol 1e-6）",
                "steps_tag": "1000 步",
            },
            {
                "file": "spring_bar.mp4",
                "name": "spring_bar  初速弹簧振荡",
                "material": "spring（边弹簧，k=youngs·l_avg=7.75e5）",
                "material_tag": "spring",
                "load": "全体初速 v=(0,0,3) m/s，无重力（纯弹簧振荡）",
                "load_tag": "初速 (0,0,3)",
                "steps": "1000 步 × dt=0.01 s（tol 1e-6）",
                "steps_tag": "1000 步",
            },
            {
                "file": "squash_sphere.mp4",
                "name": "squash_sphere  球体顶盖固定-底部下拉",
                "material": "spring（边弹簧，k=youngs·l_avg=8.20e3，sphere1K）",
                "material_tag": "spring",
                "load": "顶盖 z>0.3 固定(306 节点) + 底部 z<-0.3 区域外力 a=(0,0,-600)（318 节点），无重力",
                "load_tag": "外力 (0,0,-600)",
                "steps": "1000 步 × dt=0.01 s（tol 1e-4）",
                "steps_tag": "1000 步",
            },
        ],
    },
    "elastogen_cantilever": {
        "title": "ElastoGen cantilever 复现系列（elastogen_cantilever）",
        "series_note": "对应 rods.pdf Fig.rods-bend(b) 与 rods_twist.pdf Fig.rods-twist(a)；16×3×3 悬臂方杆（rod16x3x3.msh，144 节点/360 tet）",
        "rows": [
            {
                "file": "bend_rods_nh.mp4",
                "name": "bend_rods_nh  NH 重力弯曲",
                "material": "Neo-Hookean P1 tet FEM（E=5e7，nu=0.32）",
                "material_tag": "nh E=5e7 nu=0.32",
                "load": "重力 g=(0,0,-9.81)，左端夹持 pin-box",
                "load_tag": "重力 -9.81",
                "steps": "5000 步 × dt=0.002 s（tol 1e-3）",
                "steps_tag": "5000 步",
            },
            {
                "file": "twist_rods_nh.mp4",
                "name": "twist_rods_nh  NH 扭转 360°",
                "material": "Neo-Hookean P1 tet FEM（E=5e7，nu=0.32）",
                "material_tag": "nh E=5e7 nu=0.32",
                "load": "DBC 旋转轨迹 36°/s 绕 x 轴（10 s 恰好转 360°），无重力",
                "load_tag": "DBC 扭转 36°/s",
                "steps": "5000 步 × dt=0.002 s（tol 1e-3）",
                "steps_tag": "5000 步",
            },
            {
                "file": "bend_rods_spring.mp4",
                "name": "bend_rods_spring  弹簧对照弯曲",
                "material": "mass-spring 边弹簧（E=5e7，定性对照）",
                "material_tag": "spring E=5e7",
                "load": "重力 g=(0,0,-9.81)，左端夹持 pin-box",
                "load_tag": "重力 -9.81",
                "steps": "5000 步 × dt=0.002 s（tol 1e-3）",
                "steps_tag": "5000 步",
            },
        ],
    },
}

# 抽帧位置（帧号）与时间标注（25 fps 口径）
FRAME_SPEC = [(0, "t = 0.00 s"), (125, "t = 5.00 s"), (250, "t = 10.04 s")]


# ---------------------------------------------------------------- ffmpeg 工具
def ffprobe_frames(video):
    out = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0",
         "-show_entries", "stream=nb_frames", "-of", "csv=p=0", video],
        capture_output=True, text=True).stdout.strip()
    try:
        return int(out)
    except ValueError:
        return -1


def extract_frame(video, n, out_png):
    r = subprocess.run(
        ["ffmpeg", "-v", "error", "-y", "-i", video, "-vf", f"select=eq(n\\,{n})",
         "-vframes", "1", out_png], capture_output=True, text=True)
    return r.returncode == 0 and os.path.exists(out_png)


# ---------------------------------------------------------------- 绘制
def fit_font(draw, text, font_path_size, max_w, bold=False):
    size = font_path_size
    f = load_font(size, bold)
    while size > 12 and draw.textlength(text, font=f) > max_w:
        size -= 2
        f = load_font(size, bold)
    return f


def draw_badge(cell, lines):
    """每格左上角参数角标：半透明深色圆角块 + 两行小字。"""
    d = ImageDraw.Draw(cell, "RGBA")
    f = load_font(17)
    tw = max(d.textlength(t, font=f) for t in lines)
    bw, bh = int(tw) + 22, 24 * len(lines) + 12
    d.rectangle([8, 8, 8 + bw, 8 + bh], fill=BADGE_BG + (205,))
    d.rectangle([8, 8, 8 + bw, 8 + bh], outline=PANEL_EDGE + (255,), width=1)
    for i, t in enumerate(lines):
        d.text((19, 14 + 24 * i), t, font=f, fill=TEXT + (255,))
    return cell


def build_sheet(series, cfg, videos_dir, out_path, tmpdir):
    rows = cfg["rows"]
    ncols = len(FRAME_SPEC)
    W = MARGIN * 2 + LABEL_W + ncols * CELL_W + (ncols - 1) * GUTTER
    H = MARGIN * 2 + HEADER_H + len(rows) * (CELL_H + GUTTER) + 46
    sheet = Image.new("RGB", (W, H), BG)
    d = ImageDraw.Draw(sheet)

    # 标题
    ft = load_font(34, bold=True)
    d.text((MARGIN, MARGIN), cfg["title"], font=ft, fill=TEXT)
    fnote = load_font(19)
    d.text((MARGIN, MARGIN + 42), f"{cfg['series_note']}    生成：{date.today().isoformat()}    make_contact_sheet.py --series {series}",
           font=fnote, fill=TEXT_DIM)

    # 列头
    y0 = MARGIN + HEADER_H
    fhead = load_font(24, bold=True)
    d.text((MARGIN, y0 - 40), "实验 / 参数", font=fhead, fill=ACCENT)
    for j, (n, tlabel) in enumerate(FRAME_SPEC):
        x = MARGIN + LABEL_W + j * (CELL_W + GUTTER)
        d.rectangle([x, y0 - 46, x + CELL_W, y0 - 12], fill=PANEL)
        d.text((x + 14, y0 - 42), f"关键帧 f{n:03d}", font=fhead, fill=TEXT)
        d.text((x + CELL_W - 150, y0 - 38), tlabel, font=load_font(20), fill=TEXT_DIM)

    # 行
    for i, row in enumerate(rows):
        y = y0 + i * (CELL_H + GUTTER)
        # 行标签栏
        d.rectangle([MARGIN, y, MARGIN + LABEL_W, y + CELL_H], fill=PANEL, outline=PANEL_EDGE, width=1)
        fname = row["name"]
        f = fit_font(d, fname, 27, LABEL_W - 36, bold=True)
        d.text((MARGIN + 18, y + 16), fname, font=f, fill=TEXT)
        d.line([MARGIN + 18, y + 58, MARGIN + LABEL_W - 18, y + 58], fill=PANEL_EDGE, width=1)
        fparam = load_font(19)
        paras = [
            ("材料", row["material"]),
            ("载荷", row["load"]),
            ("步数", row["steps"]),
        ]
        yy = y + 70
        for k, v in paras:
            d.text((MARGIN + 18, yy), k, font=fparam, fill=ACCENT)
            # 参数值换行绘制（按像素宽度折行）
            words, line = v, ""
            for ch in words:
                if d.textlength(line + ch, font=fparam) > LABEL_W - 130:
                    d.text((MARGIN + 78, yy), line, font=fparam, fill=TEXT)
                    yy += 26
                    line = ch
                else:
                    line += ch
            if line:
                d.text((MARGIN + 78, yy), line, font=fparam, fill=TEXT)
                yy += 30
        d.text((MARGIN + 18, y + CELL_H - 34), row["file"], font=load_font(18), fill=TEXT_DIM)

        # 视频帧格
        video = os.path.join(videos_dir, row["file"])
        nframes = ffprobe_frames(video) if os.path.exists(video) else -1
        for j, (n, _) in enumerate(FRAME_SPEC):
            x = MARGIN + LABEL_W + j * (CELL_W + GUTTER)
            d.rectangle([x, y, x + CELL_W, y + CELL_H], outline=PANEL_EDGE, width=1)
            if nframes < 0:   # 视频缺失
                d.rectangle([x + 1, y + 1, x + CELL_W - 1, y + CELL_H - 1], fill=(46, 50, 58))
                d.text((x + 20, y + CELL_H // 2 - 12), f"缺失：{row['file']}", font=load_font(22), fill=TEXT_DIM)
                continue
            n_use, tag = n, None
            if n >= nframes:  # 帧数不齐：裁齐到末帧并标注，不留黑格
                n_use, tag = nframes - 1, f"裁齐：仅 {nframes} 帧，取末帧 f{n_use:03d}"
            png = os.path.join(tmpdir, f"{series}_{i}_{j}.png")
            if extract_frame(video, n_use, png):
                im = Image.open(png).convert("RGB").resize((CELL_W, CELL_H), Image.LANCZOS)
                cell = im.copy()
                draw_badge(cell, [
                    f"材料 {row['material_tag']}   载荷 {row['load_tag']}",
                    f"{row['steps_tag']}   帧 f{n_use:03d}" + (f"（{tag}）" if tag else ""),
                ])
                sheet.paste(cell, (x, y))
            else:
                d.rectangle([x + 1, y + 1, x + CELL_W - 1, y + CELL_H - 1], fill=(46, 50, 58))
                d.text((x + 20, y + CELL_H // 2 - 12), f"抽帧失败 f{n_use:03d}", font=load_font(22), fill=TEXT_DIM)
    return sheet


def main():
    ap = argparse.ArgumentParser(description="生成带标注的实验视频 contact sheet")
    ap.add_argument("--series", choices=sorted(MANIFEST), help="实验系列")
    ap.add_argument("--all", action="store_true", help="生成全部系列")
    ap.add_argument("--videos-dir", help="视频目录（默认 output/videos/<series>）")
    ap.add_argument("--out", help="输出 PNG（默认 output/videos/<series>/<series>_sheet.png）")
    args = ap.parse_args()

    if not args.all and not args.series:
        ap.error("需要 --series 或 --all")

    series_list = sorted(MANIFEST) if args.all else [args.series]
    tmpdir = tempfile.mkdtemp(prefix="sheet_", dir="/tmp/opencode" if os.path.isdir("/tmp/opencode") else None)
    try:
        for series in series_list:
            cfg = MANIFEST[series]
            videos_dir = args.videos_dir or os.path.join(REPO, "output", "videos", series)
            out_path = args.out or os.path.join(videos_dir, f"{series}_sheet.png")
            sheet = build_sheet(series, cfg, videos_dir, out_path, tmpdir)
            sheet.save(out_path)
            print(f"[make_contact_sheet] {out_path}  ({sheet.size[0]}x{sheet.size[1]})")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    main()
