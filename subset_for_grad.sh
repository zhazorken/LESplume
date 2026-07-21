#!/bin/bash -l
# subset_for_grad.sh — run ON CASPER. Extracts u,v,w,T,S from the time-averaged 3-D fields,
# cropped to a box around the plume, into small files for local analysis. Prints the rsync line.
#
#   ./subset_for_grad.sh                          # default box, both cases
#   X1=250 Y0=-100 Y1=100 ./subset_for_grad.sh    # tighter box (smaller files)
#   TLAST=1 ./subset_for_grad.sh                  # keep only the LAST 15-min average window
#
# Region: x = distance from ice base [m], y = along-glacier [m]; all depths (z) are kept.
set -u
module load nco 2>/dev/null || true
OUTDIR="${OUTDIR:-/glade/work/$USER/LESplume_runs}"
DEST="${DEST:-/glade/work/$USER/LESplume_sub}"; mkdir -p "$DEST"
X0=${X0:-0}; X1=${X1:-400}; Y0=${Y0:--150}; Y1=${Y1:-150}    # plume box (edit as needed)
TLAST=${TLAST:-0}                                            # 1 = keep only the last time record

for name in cg_vertical cg_overcut634; do
    src="$OUTDIR/${name}_timeavg.nc"
    [ -f "$src" ] || { echo "skip $name (no $src)"; continue; }
    out="$DEST/${name}_timeavg_plume.nc"
    tsel=""
    if [ "$TLAST" = "1" ]; then
        nt=$(ncks --trd -m -v time "$src" | grep -E 'time = .*size' | grep -oE 'size = [0-9]+' | grep -oE '[0-9]+' | tail -1)
        [ -n "$nt" ] && tsel="-d time,$((nt-1)),$((nt-1))"
    fi
    # subset every staggered horizontal dim (u on xF/yC, v on xC/yF, w on xC/yC, T,S on xC/yC)
    ncks -O -4 -L 1 -v u,v,w,T,S $tsel \
        -d xC,${X0}.,${X1}. -d xF,${X0}.,${X1}. \
        -d yC,${Y0}.,${Y1}. -d yF,${Y0}.,${Y1}. \
        "$src" "$out"
    echo "wrote $out   ($(du -h "$out" | cut -f1))"
done

echo
echo "Pull to your laptop:"
echo "  rsync -avh 'kenzhao@data-access.ucar.edu:$DEST/*_timeavg_plume.nc' ~/Desktop/Ovall26/newLES_cg/output/"
