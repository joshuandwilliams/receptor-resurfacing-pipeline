#!/usr/bin/env bash
# Compose the two panels of the RFdiffusion design-region figure from the
# ChimeraX renders in renders/.
#
# The width/height pairs are the drawn box in inches before cropping. The
# original invocation was not recorded, so these were solved back from the
# published panels' pixel dimensions and the builder's own layout constants
# (MARGIN, LABEL_PAD, COL_GAP, INPUT_SCALE, LEGEND_MIN_W, 300 dpi). They
# reproduce the published panels to within one pixel of rounding.
#
# The three design panels share one box so that the rebuilt design region 1
# length stays comparable between them; the shorter designs simply draw
# shorter after their transparent border is cropped.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p panels
python3 build_design_region_figure.py \
    --input renders/input.png 7.919 16.933 \
    --design "Design 18" "Design region 1 length: 10" renders/design_18_len10.png 9.0 18.546 \
    --design "Design 40" "Design region 1 length: 15" renders/design_40_len15.png 9.0 18.546 \
    --design "Design 48" "Design region 1 length: 20" renders/design_48_len20.png 9.0 18.546 \
    --out-a panels/rfdiff_design_regions_a.png \
    --out-b panels/rfdiff_design_regions_b.png
