#!/usr/bin/env python3
"""Check the paired article and bundled assets without running a simulation."""
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
BLOG = ROOT / 'blog'
zh = (BLOG / 'fiction.md').read_text(encoding='utf-8')
en = (BLOG / 'fiction-en.md').read_text(encoding='utf-8')

def require(condition, message):
    if not condition:
        raise SystemExit(message)

def math(source):
    return re.findall(r'\$\$[\s\S]*?\$\$|\$[^$\n]+\$', source)

require(Counter(math(zh)) == Counter(math(en)), 'Mathematical expressions differ between languages')
require(len(math(zh)) == 168, 'Unexpected mathematical expression count')
strong = lambda source: re.findall(r'\*\*(.+?)\*\*', source)
require(len(strong(zh)) == len(strong(en)) == 18, 'Emphasis spans are missing')
for a, b in zip(strong(zh)[1:-1], strong(en)[1:-1]):
    require(re.findall(r'\d+(?:\.\d+)?', a) == re.findall(r'\d+(?:\.\d+)?', b), 'Emphasized numerical values differ')
require(re.findall(r'^\*图 (\d+)：', zh, re.M) == re.findall(r'^\*Figure (\d+):', en, re.M) == [str(i) for i in range(1, 22)], 'Figure numbering differs')
images = lambda source: re.findall(r'!\[[^\]]*\]\(([^)]+)\)', source)
require(images(zh) == images(en) and len(images(zh)) == 21, 'Image references differ')
require(re.findall(r'^#{2,3} ', zh, re.M) == re.findall(r'^#{2,3} ', en, re.M), 'Heading levels differ')
for source, lang in [(zh, 'zh-CN'), (en, 'en')]:
    require(f'lang: {lang}' in source and 'draft: true' in source, 'Language or draft metadata differs')
    require('translationKey: fiction-solid-simulator' in source, 'Missing shared translation key')
    for block in source.split('\n\n')[1:]:
        if not block or block.startswith(('#', '|', '![', '$$', '```', '<!--', '*图 ', '*Figure ')):
            continue
        require(block.startswith('\u3000\u3000') if lang == 'zh-CN' else not block[0].isspace(), f'Paragraph alignment differs for {lang}')

links = 0
for file in [ROOT / 'README.md', ROOT / 'README.zh-CN.md', ROOT / 'EXPERIMENTS.md', ROOT / 'EXPERIMENTS.zh-CN.md', *BLOG.glob('*.md')]:
    source = file.read_text(encoding='utf-8')
    for href in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', source):
        if urlsplit(href).scheme or href.startswith('#'):
            continue
        require(not href.startswith('/'), f'Machine or site-specific link in {file}: {href}')
        target = file.parent / unquote(href.split('#', 1)[0])
        require(target.exists(), f'Missing link in {file.relative_to(ROOT)}: {href}')
        links += 1

manifest = BLOG / 'MANIFEST.sha256'
if manifest.exists():
    for line in manifest.read_text().splitlines():
        expected, relative = line.split('  ', 1)
        require(hashlib.sha256((BLOG / relative).read_bytes()).hexdigest() == expected, f'Hash mismatch: {relative}')

print(json.dumps({'status': 'passed', 'math_expressions_per_language': 168, 'bold_spans_per_language': 18,
                  'figures_per_language': 21, 'local_links': links, 'chinese_two_character_indent': True,
                  'english_flush_left': True, 'simulations_run': 0}, indent=2))
