# SGD-plume LES — geometric terminus, Oceananigans 0.109 (CG solver)

The newer-Oceananigans build of `../newLES`. The terminus is represented **geometrically**
(immersed ice face) with **vertical gravity** — the faithful setup — instead of the
gravity-tilt idealisation. This keeps overcut plumes ice-attached (as observed / in Ovall
et al. 2025) and removes the near-boundary pressure artifacts, because the immersed boundary
is solved with the **`ConjugateGradientPoissonSolver`** (not the approximate FFT solver).

Built on the pinned **0.109.2** environment from the sibling `IceShelfCavity` project
(`Project.toml` / `Manifest.toml` copied from there).

## What's different vs. `../newLES` (0.90.7, gravity-tilt)

| | `newLES` (0.90.7) | `newLES_cg` (0.109) |
|---|---|---|
| Terminus tilt | rotate gravity (frame tilt) | **geometry**: immersed ice face, gravity vertical |
| Floor & surface | cavity mask in a tilted box | **domain boundaries** z=0 / z=Lz (no cavity) |
| Overcut plume | tends to **detach** (adverse buoyancy) | **stays attached** (faithful) |
| Pressure (immersed) | FFT (approximate, near-wall noise) | **ConjugateGradientPoissonSolver** |
| Ambient | along ζ = z·cosθ − x·sinθ | simply **T(z)** (z = true height) |
| Quick-look | un-rotates to true x-z | already true x-z (no rotation) |

## Terminus geometry

Ice = immersed solid where `x < x_face(z) = xf_a + xf_b·z`, with β = 90 − face_angle:

- **Overcut** (`--terminus=overcut`, default): `x_face(z) = (Lz−z)·tanβ` — base sticks out to
  x = Lz·tanβ, recedes to x=0 at the surface (slopes **away** from the ocean going up).
- **Undercut** (`--terminus=undercut`): `x_face(z) = z·tanβ` — top overhangs toward the ocean.
- **Vertical** (`--face_angle=90`): no immersed ice; the ice is the west wall at x=0.

Discharge (150 m³/s, 24×6 m outlet) enters at the grounding-line corner `(x_face(0), 0)` and
rises along the ice face. Floor (z=0) and surface (z=Lz) are flat domain boundaries.

## Run it

```bash
# 0. one-time environment (needs Julia ≥ 1.10) — DIFFERENT env from newLES. There is no
#    Manifest.toml on purpose (the IceShelfCavity one was resolved on Julia 1.12 and is
#    incompatible with 1.10); your Julia resolves a fresh 0.109-compatible manifest here:
julia --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
#    (if it complains, run  Pkg.resolve()  first, or  Pkg.add("Oceananigans")  to pull 0.109.x)

# 1. CPU smoke test — compiles, steps, CG solver converges, ice geometry looks right
julia --project iceplume.jl --arch=cpu --face_angle=63.4 --stop_time=1 --simname=cgtest
python3 plot_quicklook.py cgtest      # -> output/cgtest_slices.png (immersed ice blanked)

# 2. production (GPU) — via the batch script (applies the GPU-only domain: 743 m wide ×
#    4 km long, 0.75 m fine to 375 m then stretched). Instantiate the env on the cluster first.
qsub -v CASE=vertical submit_pbs.sh
qsub -v CASE=overcut  submit_pbs.sh
```

The CLI defaults (`Ly=192`, `Lx=500`) stay small so CPU tests are cheap; the wide/long
**paper-scale domain is applied only in `submit_pbs.sh`** via flags
(`--Ly=743 --Lx=4000 --fine_x=375 --dz=0.75 --dx_max=18.6`). That's ~153 M cells → plan on an
80 GB A100 (or use `--Ly=384` for ~79 M). Because gravity is vertical here, the 4 km length
actually captures the subducting intrusion (the gravity-tilt `newLES` overcut could not).

Flags mirror `newLES` (`--discharge`, `--outlet_w/h`, `--Lz`, `--Lx`, `--fine_x`, `--dz`,
`--stop_time`, `--output_interval`, `--wall_time_limit`, `--terminus`, `--arch`, `--simname`).

## First-run CPU checks

1. **Environment**: `Pkg.instantiate()` succeeds on the 0.109.2 Manifest (separate from newLES).
2. **CG solver converges** — the progress log should step with CFL < 0.5 and no NaNs.
3. **Ice geometry** — on `cgtest_slices.png`, the immersed ice is blanked and leans the
   intended way (overcut = recedes going up); floor and surface are flat.
4. **Attachment** — the overcut plume should hug the receding ice face (the payoff vs. the
   gravity-tilt version, which detaches).

## Fidelity to Ovall et al. (2025, JGR-Oceans) LES

| spec | Ovall 2025 | here | note |
|---|---|---|---|
| Point-plume outlet | 24 m wide × 6 m tall | **24 × 6** | ✓ match (`--outlet_w/h`) |
| Injection | channel at glacier midpoint, base | at midpoint, grounding-line corner | ✓ match |
| Terminus | vertical 90°, overcut 78.7/71.6/63.4° | same (`--face_angle`, `--terminus`) | ✓ match |
| Grid resolution | 0.75 m (fine region) | **0.75 m** (`--dz`) | ✓ match |
| Fine region width | 375 m from terminus | **375 m** (`--fine_x`) | ✓ match |
| Along-fjord stretch | up to 18.6 m | **18.6 m** (`--dx_max`) | ✓ match |
| Fjord-exit outflow | above 60 m depth, 0 below, = SGD flux | same | ✓ match |
| Run length / average | 45 min / last 15 min | **45 min / 15-min window** | ✓ match |
| SGS closure | Deardorff (1980) + Ducros (1996) | AnisotropicMinimumDissipation | ~ both LES; different scheme |
| **Discharge** | 75 m³/s | **150 m³/s** | *your 2024 value* (`--discharge`) |
| **Grounding-line depth** | 172 m | **150 m** | *your 2024 value* (`--Lz`) |
| **Ambient T/S** | 2018 BPT profile | **2024 cast** | *your 2024 value* |
| Fjord width `Ly` | 743 m | 192 m (default) | cost — widen with `--Ly=743` to match |
| Fjord length `Lx` | 8,300 m | 500 m (default) | cost/near-field — `--Lx` for more far field |

Bold rows are deliberate 2024-condition updates; the last two rows are cost-driven domain
choices — pass `--Ly=743` (and a larger `--Lx`) for a paper-faithful domain if compute allows
(a point plume can feel side walls if `Ly` is too small — check the surface umbrella).

## Performance (the CG solver is slower than the old FFT solver)

The immersed CG Poisson solve is genuinely more work than `newLES`'s FFT solver — that's the
price for clean near-wall pressure. Levers, in order of impact:

- **Timestepper = QuasiAdamsBashforth2** (set by default here): ONE pressure solve per step
  instead of THREE (the RK3 default). ~3× fewer CG solves.
- **Run on GPU for anything beyond a smoke test.** The CG + immersed-boundary stack is built
  for GPUs; CPU single-thread is only for checking it compiles/steps. Don't tune CPU speed —
  use `--arch=gpu`.
- **CG tolerance / iteration cap**: `--cg_reltol` (default 1e-5; raise to 1e-4 for fewer
  iterations, lower for a stricter solve) and `--cg_maxiter` (default 30). The solver uses the
  FFT solver as a preconditioner, so it usually converges in a few iterations.
- **Cheaper CPU tests**: smaller `--Lx`, coarser `--dz`, shorter `--stop_time`; the 10 s
  face/mid-y slice cadence can be raised for long CPU runs.
- **Threads** (`julia -t auto`) *may* help on CPU, but `IceShelfCavity/NOTES.md` warns of a
  threaded-GC crash on the immersed-BC path with newer Julia — try cautiously and verify.

## Status / caveats

- **Untested in Julia here** — written against the 0.109 API using `IceShelfCavity` as the
  reference; validate on CPU first. Likely spots to adjust if something errors: the immersed
  drag `FieldBoundaryConditions(immersed=…)` signature, and `ConjugateGradientPoissonSolver`
  options — cross-check against `../IceShelfCavity/iceshelfcavity.jl` if needed.
- Immersed **drag** uses quadratic drag on all three components (ice faces). Melt is not
  included (discharge-plume study); add a 3-equation melt BC later if wanted (the pattern is
  in `IceShelfCavity/scripts/melt_parameterization.jl`).
- Ambient fit, `ambient_profile.jl`, and `verify_setup.py` (BPT) are shared with `newLES`.
