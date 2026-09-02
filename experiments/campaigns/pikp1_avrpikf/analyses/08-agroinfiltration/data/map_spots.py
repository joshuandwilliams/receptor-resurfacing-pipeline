"""Map the slide-2 annotation boxes onto full-resolution image pixels."""
import json
import re
import zipfile
from pathlib import Path

PPTX = ("/Users/jowillia/Documents/GitHub/receptor-resurfacing-pipeline/experiments/"
        "campaigns/pikp1_avrpikf/runs/crystal_full_test_contig/"
        "AI_Pikp1_HMA_Prelim_2_JoshW.pptx")
OUT = Path(__file__).parent / "spots.json"
IMG_W, IMG_H = 5568, 3712
EMU = 914400.0

z = zipfile.ZipFile(PPTX)
x = z.read("ppt/slides/slide2.xml").decode("utf8")

pic = re.search(r"<p:pic>.*?</p:pic>", x, re.S).group(0)
sr = re.search(r'<a:srcRect l="(\d+)" t="(\d+)" r="(\d+)" b="(\d+)"/>', pic)
l, t, r, b = (int(v) / 100000.0 for v in sr.groups())
po = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"/><a:ext cx="(\d+)" cy="(\d+)"', pic)
px, py, pw, ph = (int(v) / EMU for v in po.groups())

x0, x1 = l * IMG_W, (1 - r) * IMG_W
y0, y1 = t * IMG_H, (1 - b) * IMG_H

# Walk the shapes in document order; each run of boxes ends with its legend swatch
# (at x ~9.85) followed by the text box naming that group.
groups, cur = [], []
for m in re.finditer(r"<p:sp>.*?</p:sp>", x, re.S):
    s = m.group(0)
    off = re.search(r'<a:off x="(-?\d+)" y="(-?\d+)"/><a:ext cx="(\d+)" cy="(\d+)"', s)
    txt = "".join(re.findall(r"<a:t>(.*?)</a:t>", s)).strip()
    if not off:
        continue
    ox, oy, ow, oh = (int(v) / EMU for v in off.groups())
    if txt:
        if cur:
            groups.append((txt, cur))
            cur = []
        continue
    if ox > 9.5:          # legend swatch, not a leaf annotation
        continue
    cx, cy = ox + ow / 2, oy + oh / 2
    fx, fy = (cx - px) / pw, (cy - py) / ph
    cur.append((round(x0 + fx * (x1 - x0)), round(y0 + fy * (y1 - y0))))

spots = []
for label, pts in groups:
    label = re.sub(r"\s+", " ", label)
    for p in pts:
        spots.append({"label": label, "x": p[0], "y": p[1]})

OUT.write_text(json.dumps(spots, indent=1))
for label, pts in groups:
    print(f"{len(pts):2d}  {re.sub(r'\\s+', ' ', label)[:70]}")
