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
    ap.add_argument("--xmax", type=float, default=400.0, help="near-field x-limit [m] for midy (0=full)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    fn = os.path.join(a.dir, f"{a.prefix}_{a.slice}.nc")
    if not os.path.exists(fn): sys.exit(f"not found: {fn}")
    ds = xr.open_dataset(fn, decode_timedelta=False)
    horiz = "x" if a.slice == "midy" else "y"
    xlab = "distance from ice base (m)" if a.slice == "midy" else "along-glacier y (m)"
    xf_a, xf_b, face_x = attr(ds, "xf_a"), attr(ds, "xf_b"), attr(ds, "face_x_m", 0.0)

    # Load + ice-mask EACH variable with ITS OWN coords (w is on z-faces, T/S on z-centers, so
    # they differ by one point in z — a single shared mask would mis-broadcast).
    series, coords = {}, {}
    for v in ("w", "T", "S"):
        if v not in ds: continue
        arr, Hv, Zv = load_series(ds, v, horiz)
        if a.slice == "midy":
            icev = Hv[None, :] < (xf_a + xf_b * Zv[:, None])
        else:
            icev = np.broadcast_to((face_x < (xf_a + xf_b * Zv))[:, None], (Zv.size, Hv.size))
        series[v] = np.where(icev[None] | ~np.isfinite(arr), np.nan, arr)
        coords[v] = (Hv, Zv)
    if not series: sys.exit("no w/T/S in file")
    nt = next(iter(series.values())).shape[0]
    times = _coord(ds, "time")[:nt] if "time" in ds.variables else np.arange(nt)

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
        Hv, Zv = coords[v]; lo, hi, cmap = rng[v]
        pc = ax.pcolormesh(Hv, Zv, series[v][0], cmap=cmap, vmin=lo, vmax=hi, shading="auto")
        fig.colorbar(pc, ax=ax, label={"w": "w (m/s)", "T": "T (°C)", "S": "S (g/kg)"}[v])
        ax.set_xlabel(xlab); ax.set_ylabel("height above grounding line (m)"); ax.set_title(v)
        if a.slice == "midy" and a.xmax: ax.set_xlim(0, a.xmax)   # near-field zoom
        meshes.append(pc)
    title = fig.suptitle("")
    fig.tight_layout(rect=[0, 0, 1, 0.96])       # leave room for the suptitle; de-crowd panels

    def update(i):
        for pc, v in zip(meshes, vars_present):
            pc.set_array(series[v][i].ravel())      # shading='auto'→'nearest': C is nZ×nH
        title.set_text(f"{a.prefix} — {a.slice}   t = {times[i]:.0f} s   (frame {i+1}/{nt})")
        return meshes

    ani = manim.FuncAnimation(fig, update, frames=nt, blit=False)
    base = a.out and os.path.splitext(a.out)[0] or os.path.join(a.dir, f"{a.prefix}_{a.slice}")
    if manim.FFMpegWriter.isAvailable():
        # H.264 + yuv420p + even pixel dims = plays in QuickTime/Preview/browsers
        out = base + ".mp4"
        w = manim.FFMpegWriter(fps=a.fps, codec="libx264",
                               extra_args=["-pix_fmt", "yuv420p",
                                           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2"])
        ani.save(out, writer=w, dpi=120)
    else:
        out = base + ".gif"
        print("  (ffmpeg not found — writing a .gif instead; bring THIS file back)")
        ani.save(out, writer=manim.PillowWriter(fps=a.fps), dpi=90)
    plt.close(fig); ds.close()
    sz = os.path.getsize(out) / 1e6
    print(f"  wrote {out}  ({nt} frames, {sz:.1f} MB)")
    if sz < 0.01: print("  WARNING: output is tiny — the writer likely failed; check ffmpeg/Pillow.")

if __name__ == "__main__":
    main()
