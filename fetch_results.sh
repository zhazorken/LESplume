#!/bin/bash
# fetch_results.sh — run ON YOUR LAPTOP to pull the LIGHT results from Casper into ./output:
# the 2-D slices, the time-average, and any PNGs/movies postprocess.sh made. Leaves the heavy
# 3-D *_fields.nc on GLADE (use `ncks` there to subset a slab if you need volume data locally).
#
#   ./fetch_results.sh
#   REMOTE=kenzhao@data-access.ucar.edu:/glade/work/kenzhao/LESplume_runs ./fetch_results.sh
#
# For very large / many-file transfers, prefer Globus over rsync.
set -u
REMOTE="${REMOTE:-kenzhao@data-access.ucar.edu:/glade/work/kenzhao/LESplume_runs}"
DEST="${DEST:-output}"
mkdir -p "$DEST"

rsync -avh --progress \
  --include='*_midy.nc' --include='*_face.nc' --include='*_timeavg.nc' \
  --include='*.png' --include='*.mp4' --include='*.gif' \
  --exclude='*' \
  "$REMOTE/" "$DEST/"

echo "Pulled light outputs to $DEST/. Re-plot locally with:  python3 plot_quicklook.py <prefix> --dir $DEST"
