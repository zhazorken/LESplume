#!/usr/bin/env python3
"""Regenerate the GPU-safe ambient-profile fit embedded in ambient_profile.jl.

Reads AmbientProfile_15JUL2024_4modeling.mat (fields amb.CT, amb.SA, amb.depth), converts
to the model coordinate z = Lz - depth (height above the grounding line), holds properties
constant below the moraine crest (depth 125 m => z = 25 m), and fits an endpoint-weighted
degree-8 polynomial in the normalised height zeta = clamp((z-z_crest)/(Lz-z_crest),0,1) to
the variable region. Prints the coefficients (highest power first) to paste into
ambient_profile.jl, and writes ambient_2024.csv + ambient_fit.png.

Usage:  python3 scripts/fit_ambient.py /path/to/AmbientProfile_15JUL2024_4modeling.mat
"""
import sys, csv
import numpy as np
import scipy.io as sio

MAT = sys.argv[1] if len(sys.argv) > 1 else "../AmbientProfile_15JUL2024_4modeling.mat"
LZ = 150.0            # grounding-line depth [m]
CREST_DEPTH = 125.0   # lowest point on the moraine crest [m]
DEG = 8

amb = sio.loadmat(MAT)["amb"]
CT = np.array(amb["CT"][0, 0]).flatten()
SA = np.array(amb["SA"][0, 0]).flatten()
depth = np.array(amb["depth"][0, 0]).flatten()

sel = (depth >= 0) & (depth <= LZ)
z = (LZ - depth)[sel]
T = CT[sel]; S = SA[sel]
order = np.argsort(z); z, T, S = z[order], T[order], S[order]
z_crest = LZ - CREST_DEPTH

zeta = lambda zv: np.clip((zv - z_crest) / (LZ - z_crest), 0, 1)
mz = z >= z_crest
zt = zeta(z[mz])
w = np.ones_like(zt); w[(zt < 0.06) | (zt > 0.94)] = 8.0   # weight the surface & crest ends

def wpolyfit(x, y, deg, w):
    V = np.vander(x, deg + 1); W = np.sqrt(w)[:, None]
    c, *_ = np.linalg.lstsq(V * W, y * np.sqrt(w), rcond=None)
    return c

cT = wpolyfit(zt, T[mz], DEG, w)
cS = wpolyfit(zt, S[mz], DEG, w)
Tfit = np.polyval(cT, zeta(z)); Sfit = np.polyval(cS, zeta(z))
print(f"T RMSE {np.sqrt(np.mean((Tfit-T)**2)):.4f}  max {np.max(np.abs(Tfit-T)):.4f}")
print(f"S RMSE {np.sqrt(np.mean((Sfit-S)**2)):.4f}  max {np.max(np.abs(Sfit-S)):.4f}")
print("CT_COEFFS =", ", ".join(f"{v:.8g}" for v in cT))
print("CS_COEFFS =", ", ".join(f"{v:.8g}" for v in cS))

with open("ambient_2024.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["z_model_m", "depth_m", "CT_data", "SA_data", "T_fit", "S_fit"])
    for i in range(len(z) - 1, -1, -1):
        wr.writerow([f"{z[i]:.3f}", f"{LZ-z[i]:.3f}", f"{T[i]:.4f}", f"{S[i]:.4f}",
                     f"{Tfit[i]:.4f}", f"{Sfit[i]:.4f}"])

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    ax[0].plot(T, z, ".", ms=2, color="tab:red", label="2024 data")
    ax[0].plot(Tfit, z, "-k", lw=2, label="fit")
    ax[0].axhline(z_crest, ls="--", c="gray"); ax[0].set_xlabel("CT (°C)")
    ax[0].set_ylabel("z above grounding line (m)"); ax[0].set_title("Temperature"); ax[0].legend()
    ax[1].plot(S, z, ".", ms=2, color="tab:blue", label="2024 data")
    ax[1].plot(Sfit, z, "-k", lw=2, label="fit")
    ax[1].axhline(z_crest, ls="--", c="gray", label="moraine crest")
    ax[1].set_xlabel("SA (g/kg)"); ax[1].set_title("Salinity"); ax[1].legend()
    fig.suptitle("2024 ambient profile (constant below moraine crest) — data vs GPU-safe fit")
    fig.tight_layout(); fig.savefig("ambient_fit.png", dpi=110)
except Exception as e:
    print("plot skipped:", e)
