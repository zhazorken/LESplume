#!/usr/bin/env python3
"""Animate the SGD-plume slice time-series (geometric / vertical-gravity output).

The `*_midy.nc` (x-z) and `*_face.nc` (y-z) slices are written every ~10 s, so they make a
good movie of the plume spin-up. This renders a 3-panel (w, T, S) animation to mp4 (ffmpeg)
or gif (fallback), with the immersed ice blanked and color scales fixed across frames.

Usage:
  python3 make_movie.py cg_overcut634                    # midy slice, both->mp4
  python3 make_movie.py cg_overcut634 --slice face --dir /glade/work/$USER/LESplume_runs
  python3 make_movie.py cg_vertical --fps 12 --out vert.mp4
"""
import argparse, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manim
try:
    import xarray as xr
except ImportError:
    sys.exit("Please: pip install xarray netCDF4")

def _ax(da, l):
    for d in da.dims:
        if d.lower().startswith(l): return d
    return None
def _coord(ds, d):
    return np.asarray(ds[d].values) if d in ds.variables else np.arange(ds.sizes[d])
def attr(ds, k, d=0.0):
    return float(np.ravel(ds.attrs[k])[0]) if k in ds.attrs else d

def load_series(ds, v, horiz):
    """Return (data[t,z,h], Hcoord, Zcoord) for a slice variable over all times."""
    da = ds[v].squeeze()
    hd, zd = _ax(da, horiz), _ax(da, "z")
    td = "time" if "time" in da.dims else None
    if td is None:
        arr = da.transpose(zd, hd).values[None]           # single frame
    else:
        arr = da.transpose(td, zd, hd).values
    return arr, _coord(ds, hd), _coord(ds, zd)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix")
    ap.add_argument("--dir", default="output")
    ap.add_argument("--slice", default="midy", choices=["midy", "face"])
    ap.add_argument("--fps", type=int, default=10)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    fn = os.path.join(a.dir, f"{a.prefix}_{a.slice}.nc")
    if not os.path.exists(fn): sys.exit(f"not found: {fn}")
    ds = xr.open_dataset(fn)
    horiz = "x" if a.slice == "midy" else "y"
    xlab = "distance from ice base (m)" if a.slice == "midy" else "along-glacier y (m)"
    xf_a, xf_b, face_x = attr(ds, "xf_a"), attr(ds, "xf_b"), attr(ds, "face_x_m", 0.0)

    series, H, Z = {}, None, None
    for v in ("w", "T", "S"):
        if v in ds:
            series[v], H, Z = load_series(ds, v, horiz)
    if not series: sys.exit("no w/T/S in file")
    nt = next(iter(series.values())).shape[0]
    times = _coord(ds, "time")[:nt] if "time" in ds.variables else np.arange(nt)

    # immersed-ice mask (same every frame): solid where x < x_face(z)
    if a.slice == "midy":
        ice = H[None, :] < (xf_a + xf_b * Z[:, None])
    else:
        ice = np.broadcast_to((face_x < (xf_a + xf_b * Z))[:, None], (Z.size, H.size))
    for v in series:
        series[v] = np.where(ice[None] | ~np.isfinite(series[v]), np.nan, series[v])

    # fixed color scales over the whole series (robust percentiles)
    rng = {}
    for v, arr in series.items():
        if v == "w":
            m = np.nanpercentile(np.abs(arr), 99) or 1e-6; rng[v] = (-m, m, "RdBu_r")
        else:
            lo, hi = np.nanpercentile(arr, [1, 99]); rng[v] = (lo, hi, "inferno" if v == "T" else "viridis")

    vars_present = [v for v in ("w", "T", "S") if v in series]
    fig, axs = plt.subplots(1, len(vars_present), figsize=(5.2 * len(vars_present), 5), squeeze=False)
    axs = axs[0]
    meshes = []
    for ax, v in zip(axs, vars_present):
        lo, hi, cmap = rng[v]
        pc = ax.pcolormesh(H, Z, series[v][0], cmap=cmap, vmin=lo, vmax=hi, shading="auto")
        fig.colorbar(pc, ax=ax, label={"w": "w (m/s)", "T": "T (°C)", "S": "S (g/kg)"}[v])
        ax.set_xlabel(xlab); ax.set_ylabel("height above grounding line (m)"); ax.set_title(v)
        meshes.append(pc)
    title = fig.suptitle("")
    fig.tight_layout(rect=[0, 0, 1, 0.96])       # leave room for the suptitle; de-crowd panels

    def update(i):
        for pc, v in zip(meshes, vars_present):
            pc.set_array(series[v][i].ravel())      # shading='auto'→'nearest': C is nZ×nH
        title.set_text(f"{a.prefix} — {a.slice}   t = {times[i]:.0f} s   (frame {i+1}/{nt})")
        return meshes

    ani = manim.FuncAnimation(fig, update, frames=nt, blit=False)
    out = a.out or os.path.join(a.dir, f"{a.prefix}_{a.slice}.mp4")
    try:
        ani.save(out, writer=manim.FFMpegWriter(fps=a.fps, bitrate=3000))
    except Exception as e:
        out = os.path.splitext(out)[0] + ".gif"
        print(f"  (ffmpeg unavailable [{e}] — writing gif instead)")
        ani.save(out, writer=manim.PillowWriter(fps=a.fps))
    plt.close(fig)
    print(f"  wrote {out}  ({nt} frames)")
    ds.close()

if __name__ == "__main__":
    main()
