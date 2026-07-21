#!/bin/bash -l
# subset_for_grad.sh — run ON CASPER. Extract u,v,w,T,S from the time-averaged 3-D fields,
# cropped to a box around the plume, into small compressed files for local analysis.
# Uses xarray (not ncks) so it works regardless of Oceananigans' NetCDF dimension names.
#
#   ./subset_for_grad.sh                          # default box, both cases
#   X1=250 Y0=-100 Y1=100 ./subset_for_grad.sh    # tighter box -> smaller files
#   TLAST=1 ./subset_for_grad.sh                  # keep only the LAST 15-min average window
#
# Region: x = distance from ice base [m], y = along-glacier [m]; all depths (z) kept.
set -u
module load conda 2>/dev/null && conda activate npl 2>/dev/null || echo "(assuming python w/ xarray on PATH)"
PY="${PY:-python}"
export OUTDIR="${OUTDIR:-/glade/work/$USER/LESplume_runs}"
export DEST="${DEST:-/glade/work/$USER/LESplume_sub}"; mkdir -p "$DEST"
export X0=${X0:-0} X1=${X1:-400} Y0=${Y0:--150} Y1=${Y1:-150} TLAST=${TLAST:-0}

$PY - <<'PYEOF'
import os, xarray as xr
X0,X1,Y0,Y1 = (float(os.environ[k]) for k in ("X0","X1","Y0","Y1"))
TLAST = os.environ["TLAST"] == "1"; OUTDIR, DEST = os.environ["OUTDIR"], os.environ["DEST"]
for name in ("cg_vertical", "cg_overcut634"):
    src = os.path.join(OUTDIR, f"{name}_timeavg.nc")
    if not os.path.exists(src):
        print("skip", name, "(no", src + ")"); continue
    ds = xr.open_dataset(src, decode_timedelta=False)[["u", "v", "w", "T", "S"]]
    # crop every x-dim to [X0,X1] and every y-dim to [Y0,Y1], whatever they're named
    sel = {}
    for d in ds.dims:
        if d in ds.variables or d in ds.coords:
            if d[0].lower() == "x": sel[d] = slice(X0, X1)
            elif d[0].lower() == "y": sel[d] = slice(Y0, Y1)
    ds = ds.sel(**sel)
    if TLAST and "time" in ds.dims: ds = ds.isel(time=[-1])
    out = os.path.join(DEST, f"{name}_timeavg_plume.nc")
    enc = {v: {"zlib": True, "complevel": 1} for v in ds.data_vars}
    ds.to_netcdf(out, encoding=enc)
    print(f"wrote {out}  {os.path.getsize(out)/1e6:.0f} MB  dims={dict(ds.sizes)}")
    ds.close()
PYEOF

echo
echo "Pull to your laptop:"
echo "  rsync -avh 'kenzhao@data-access.ucar.edu:$DEST/*_timeavg_plume.nc' ~/Desktop/Ovall26/newLES_cg/output/"
