#!/usr/bin/env python3
"""Quick-look plots for the GEOMETRIC (immersed-ice, vertical-gravity) SGD-plume LES.

Because gravity is vertical, the output is already in true x-z — no rotation. The immersed
ice (x < x_face(z) = xf_a + xf_b·z, read from NetCDF metadata) is blanked. Produces:
  <prefix>_slices.png   mid-fjord (x-z) and near-grounding-line (y-z) snapshots of w, T, S
  <prefix>_profiles.png upward volume flux Q(z), flux-weighted T(z)/S(z), mean/max w(z),
                        with the neutral-buoyancy level

Usage:  python3 plot_quicklook.py [prefix] [--dir output] [--time -1]
"""
import argparse, glob, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
try:
    import xarray as xr
except ImportError:
    sys.exit("Please: pip install xarray netCDF4")

ALPHA, BETA, RHO0 = 3.87e-5, 7.86e-4, 1027.0

def _ax(da, l):
    for d in da.dims:
        if d.lower().startswith(l): return d
    return None
def _coord(ds, d):
    return np.asarray(ds[d].values) if d in ds.variables else np.arange(ds.sizes[d])
def _w(c):
    c = np.asarray(c, float)
    if c.size < 2: return np.array([1.0])
    e = np.empty(c.size+1); e[1:-1] = 0.5*(c[:-1]+c[1:]); e[0]=c[0]-(c[1]-c[0])/2; e[-1]=c[-1]+(c[-1]-c[-2])/2
    return np.diff(e)
def attr(ds, k, d=0.0):
    return float(np.ravel(ds.attrs[k])[0]) if k in ds.attrs else d

def _load2d(ds, v, tidx, horiz):
    if v not in ds: return None
    da = ds[v].isel(time=tidx) if "time" in ds[v].dims else ds[v]
    da = da.squeeze()
    hd, zd = _ax(da, horiz), _ax(da, "z")
    if hd is None or zd is None: return None
    return da.transpose(zd, hd).values, _coord(ds, hd), _coord(ds, zd)

def plot_slices(prefix, outdir, tidx, xmax=None):
    fig, axs = plt.subplots(2, 3, figsize=(15, 8))
    rows = [("midy", "x", "distance from ice base (m)"),
            ("face", "y", "along-glacier y (m)")]
    for r, (tag, horiz, xlab) in enumerate(rows):
        fn = os.path.join(outdir, f"{prefix}_{tag}.nc")
        if not os.path.exists(fn):
            for c in range(3): axs[r, c].set_visible(False)
            print(f"  (skip {tag}: {fn} not found)"); continue
        ds = xr.open_dataset(fn, decode_timedelta=False)
        xf_a, xf_b = attr(ds, "xf_a"), attr(ds, "xf_b")
        face_x = attr(ds, "face_x_m", 0.0)
        for c, (v, cmap) in enumerate([("w", "RdBu_r"), ("T", "inferno"), ("S", "viridis")]):
            ax = axs[r, c]
            got = _load2d(ds, v, tidx, horiz)
            if got is None: ax.set_visible(False); continue
            data, H, Z = got
            # blank the immersed ice (x < x_face(z)); face slice is at fixed x=face_x
            if tag == "midy":
                ice = H[None, :] < (xf_a + xf_b * Z[:, None])
            else:
                ice = np.broadcast_to((face_x < (xf_a + xf_b * Z))[:, None], data.shape)
            data = np.where(ice | ~np.isfinite(data), np.nan, data)
            if v == "w":
                m = np.nanpercentile(np.abs(data), 99) or 1e-6
                pc = ax.pcolormesh(H, Z, data, cmap=cmap, vmin=-m, vmax=m, shading="auto")
            else:
                lo, hi = np.nanpercentile(data, [1, 99])
                pc = ax.pcolormesh(H, Z, data, cmap=cmap, vmin=lo, vmax=hi, shading="auto")
            fig.colorbar(pc, ax=ax, label={"w":"w (m/s)","T":"T (°C)","S":"S (g/kg)"}[v])
            ax.set_xlabel(xlab); ax.set_ylabel("height above grounding line (m)")
            ax.set_title(f"{tag}: {v}")
            if tag == "midy" and xmax: ax.set_xlim(0, xmax)   # near-field zoom
        ds.close()
    fig.suptitle(f"{prefix} — snapshot slices (time index {tidx})  [true x-z, immersed ice]")
    fig.tight_layout()
    out = os.path.join(outdir, f"{prefix}_slices.png"); fig.savefig(out, dpi=120); plt.close(fig)
    print("  wrote", out)

def plot_profiles(prefix, outdir):
    tavg = os.path.join(outdir, f"{prefix}_timeavg.nc")
    fields = os.path.join(outdir, f"{prefix}_fields.nc")
    ds = None
    if os.path.exists(tavg):
        d = xr.open_dataset(tavg, decode_timedelta=False)
        if "w" in d and d["w"].size > 0 and float(np.abs(d["w"]).max()) > 1e-6: ds = d
        else: d.close(); print("  (timeavg empty — using instantaneous fields)")
    if ds is None and os.path.exists(fields): ds = xr.open_dataset(fields, decode_timedelta=False)
    if ds is None: print("  (skip profiles: no output)"); return
    if "time" in ds.coords or "time" in ds.variables:
        tval = float(np.ravel(ds["time"].values)[-1])
        print(f"  profiles from snapshot t = {tval:.0f} s")
        if tval < 30:
            print("  WARNING: that's essentially the INITIAL condition — the run is too short "
                  "for the 3D fields/time-average cadence. The slices (10 s) show the real "
                  "state; lengthen the run or lower --output_interval for meaningful profiles.")

    def arr(v):
        v = v.isel(time=-1) if "time" in v.dims else v
        return v.transpose(_ax(v, "x"), _ax(v, "y"), _ax(v, "z")).values
    W, Tv, Sv = arr(ds["w"]), arr(ds["T"]), arr(ds["S"])
    nx, ny, nz = (min(a.shape[i] for a in (W, Tv, Sv)) for i in range(3))
    W, Tv, Sv = (a[:nx, :ny, :nz] for a in (W, Tv, Sv))
    x = _coord(ds, _ax(ds["w"], "x"))[:nx]; z = _coord(ds, _ax(ds["w"], "z"))[:nz]
    area = np.outer(_w(x), _w(_coord(ds, _ax(ds["w"], "y"))[:ny]))
    xf_a, xf_b = attr(ds, "xf_a"), attr(ds, "xf_b")

    # exclude immersed-ice cells (x < x_face(z)) and any NaNs
    ice = (x[:, None] < (xf_a + xf_b * z[None, :]))[:, None, :]     # (x,1,z)
    solid = np.broadcast_to(ice, W.shape) | ~np.isfinite(W)
    Wm = np.where(solid, np.nan, W); Tm = np.where(solid, np.nan, Tv); Sm = np.where(solid, np.nan, Sv)

    A = np.broadcast_to(area[:, :, None], W.shape)
    far = slice(int(nx * 0.66), None)
    with np.errstate(invalid="ignore", divide="ignore"):
        Tamb = np.nanmean(Tm[far], axis=(0, 1)); Samb = np.nanmean(Sm[far], axis=(0, 1))
        # robust per-level ambient S (saltiest ~90th pct = unmixed water); PLUME = fresher by dS.
        # Integrating flux only over plume cells stops ambient turbulence from inflating Q (the
        # surfaced vertical plume otherwise makes the whole-domain flux meaningless).
        Sref = np.nanpercentile(np.where(np.isfinite(Sm), Sm, np.nan), 90, axis=(0, 1))
    dS = 0.1                                   # g/kg freshness threshold defining "plume water"
    plume = np.isfinite(Wm) & (Sm < (Sref[None, None, :] - dS))
    Wp = np.where(plume, np.nan_to_num(Wm), 0.0)
    up = np.where(Wp > 0, Wp, 0.0)             # upward-moving plume water only
    Ap = np.where(plume, A, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        Q = np.nansum(A * up, axis=(0, 1))                                      # plume volume flux
        Tbar = np.nansum(A * up * np.nan_to_num(Tm), axis=(0, 1)) / np.where(Q > 0, Q, np.nan)
        Sbar = np.nansum(A * up * np.nan_to_num(Sm), axis=(0, 1)) / np.where(Q > 0, Q, np.nan)
        Apz = Ap.sum((0, 1))
        wmean = np.nansum(Ap * Wp, axis=(0, 1)) / np.where(Apz > 0, Apz, np.nan)  # plume-mean w
        wmax = np.nanmax(np.where(plume, Wm, -np.inf), axis=(0, 1)); wmax[~np.isfinite(wmax)] = np.nan

    rho_p = RHO0 * (1 - ALPHA * Tbar + BETA * Sbar)
    rho_a = RHO0 * (1 - ALPHA * Tamb + BETA * Samb)
    nb = None
    idx = np.where(np.diff(np.sign(rho_p - rho_a)) > 0)[0]
    if idx.size: nb = z[idx[0] + 1]

    fig, ax = plt.subplots(1, 4, figsize=(17, 6), sharey=True)
    ax[0].plot(Q, z); ax[0].set_xlabel("plume volume flux Q(z) [m³/s]"); ax[0].set_ylabel("height above grounding line (m)")
    ax[1].plot(Tbar, z, label="plume (flux-wtd)"); ax[1].plot(Tamb, z, "--", c="gray", label="ambient"); ax[1].set_xlabel("T (°C)"); ax[1].legend()
    ax[2].plot(Sbar, z, label="plume (flux-wtd)"); ax[2].plot(Samb, z, "--", c="gray", label="ambient"); ax[2].set_xlabel("S (g/kg)"); ax[2].legend()
    ax[3].plot(wmean, z, label="plume-mean w"); ax[3].plot(wmax, z, label="max w (plume)"); ax[3].set_xlabel("w (m/s)"); ax[3].legend()
    for a in ax:
        a.grid(alpha=0.3)
        if nb is not None: a.axhline(nb, ls=":", c="crimson")
    ttl = f"{prefix} — plume profiles" + (f"  (neutral buoyancy z≈{nb:.0f} m)" if nb else "")
    fig.suptitle(ttl); fig.tight_layout()
    out = os.path.join(outdir, f"{prefix}_profiles.png"); fig.savefig(out, dpi=120); plt.close(fig)
    print("  wrote", out); ds.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("prefix", nargs="?", default=None)
    ap.add_argument("--dir", default="output")
    ap.add_argument("--time", type=int, default=-1)
    ap.add_argument("--xmax", type=float, default=400.0,
                    help="near-field x-limit [m] for the midy panels (0/None = full domain)")
    a = ap.parse_args()
    prefix = a.prefix
    if prefix is None:
        cands = sorted(glob.glob(os.path.join(a.dir, "*_midy.nc")) + glob.glob(os.path.join(a.dir, "*_fields.nc")))
        if not cands: sys.exit(f"No output in {a.dir}/. Pass a prefix.")
        prefix = os.path.basename(cands[0]).rsplit("_", 1)[0]
    os.makedirs(a.dir, exist_ok=True)
    print(f"Plotting run '{prefix}' from {a.dir}/")
    plot_slices(prefix, a.dir, a.time, xmax=(a.xmax or None))
    plot_profiles(prefix, a.dir)

if __name__ == "__main__":
    main()
