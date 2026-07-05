#!/bin/bash -l
# postprocess.sh — run ON CASPER (login node or a small job) to make quick-look PNGs + movies
# from the run outputs, WHERE THE DATA LIVES (so you never move the 7 GB 3-D files just to plot).
#
#   ./postprocess.sh                       # both cases, from $OUTDIR
#   OUTDIR=/glade/derecho/scratch/$USER/LESplume ./postprocess.sh
#   ./postprocess.sh cg_overcut634         # a single run
#
# Needs Python with xarray/netCDF4/matplotlib (+ffmpeg for mp4). NCAR provides these via conda:
set -u
cd "$(dirname "$0")"
OUTDIR="${OUTDIR:-/glade/work/$USER/LESplume_runs}"

module load conda 2>/dev/null && conda activate npl 2>/dev/null || \
    echo "(couldn't load conda/npl — assuming python w/ xarray,matplotlib is already on PATH)"
PY="${PY:-python}"

RUNS=("$@"); [ ${#RUNS[@]} -eq 0 ] && RUNS=(cg_vertical cg_overcut634)
for name in "${RUNS[@]}"; do
    echo "=== $name ==="
    $PY plot_quicklook.py "$name" --dir "$OUTDIR" || echo "  (quick-look skipped)"
    $PY make_movie.py    "$name" --dir "$OUTDIR" --slice midy || echo "  (midy movie skipped)"
    $PY make_movie.py    "$name" --dir "$OUTDIR" --slice face || echo "  (face movie skipped)"
done
echo "Done. PNGs + movies are in $OUTDIR (fetch them with fetch_results.sh from your laptop)."
