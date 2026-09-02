#!/usr/bin/env bash
# Rebuild the two-panel design-region MSA figure for the thesis.
#
# The pipeline's own run_alignment.py (cached alongside the .aln files in data/)
# renders region 1 pale and region 2 saturated, and has no SNK-EKE row.  Its
# output still sits in the run folder on the HPC and is deliberately not
# touched.  The thesis figure is the re-rendered version instead, built here.
#
# Chain:
#   1. make_msa_plots.py            the alignments        -> plots/region{1,2}_msa.png
#                                                        plots/colour_key.png
#   2. assemble_figure.py       panels + key          -> design_region_msa.svg
#   3. inkscape                 the svg               -> design_region_msa.png
#   4. copy svg + png into every chapter figure directory listed below
#
# Step 2 uses the thesis's shared assembler, which adds the A/B letters and the
# A4 layout.  Step 3 exports the drawing bounding box at 300 dpi.  That box is
# the 155mm content width every thesis figure shares, which is where the common
# 1831 px width comes from.  Exporting the bbox-spacer rect instead loses the
# first panel's A label, because the letter is drawn above the rect's top edge.
#
# Usage:  ./rebuild_thesis_figure.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THESIS="${THESIS:-$HOME/Documents/thesis}"
ASSEMBLER="$THESIS/chapters/05-nextflow-pipeline-max-effectors/scripts/assemble_figure.py"

# The 0X mirror chapter is gone, so there is one destination now.
FIG_DIRS=(
    "$THESIS/chapters/05-nextflow-pipeline-max-effectors/figures/thesis-figures"
)

STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
SVG="$STAGE/design_region_msa.svg"
PNG="$STAGE/design_region_msa.png"

[[ -f "$ASSEMBLER" ]] || { echo "assembler not found: $ASSEMBLER" >&2; exit 1; }
# Inkscape ships as an app bundle on macOS and is normally not on PATH.
INKSCAPE="${INKSCAPE:-$(command -v inkscape || true)}"
[[ -n "$INKSCAPE" ]] || INKSCAPE=/Applications/Inkscape.app/Contents/MacOS/inkscape
[[ -x "$INKSCAPE" ]] || { echo "inkscape not found (set INKSCAPE=...)" >&2; exit 1; }

echo "1/4  panels"
python3 "$HERE/make_msa_plots.py"

echo "2/4  assemble"
# The key goes in as an --extra so it gets no panel letter and a tighter gap.
python3 "$ASSEMBLER" "$SVG" \
    --panel A "$HERE/plots/region1_msa.png" \
    --panel B "$HERE/plots/region2_msa.png" \
    --extra "$HERE/plots/colour_key.png"

echo "3/4  export png"
"$INKSCAPE" "$SVG" \
    --export-area-drawing \
    --export-dpi=300 \
    --export-type=png \
    --export-filename="$PNG" >/dev/null

python3 -c "
from PIL import Image
import sys
w, h = Image.open('$PNG').size
print(f'     {w} x {h} px')
if w != 1831:
    sys.exit(f'expected 1831 px wide, got {w}')
"

echo "4/4  install"
for d in "${FIG_DIRS[@]}"; do
    if [[ -d "$d" ]]; then
        cp "$SVG" "$PNG" "$d/"
        echo "     $d"
    else
        echo "     skipped (missing): $d"
    fi
done
echo "done"
