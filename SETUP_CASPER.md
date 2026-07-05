# Getting these runs onto Casper (NCAR) via GitHub

The code is self-contained (the ambient profile is baked into `ambient_profile.jl`), so a small
GitHub repo is the easiest transfer. No Manifest or data files are committed — the environment
resolves fresh on the cluster.

## 1. Create the repo and push (from your laptop, in this folder)

Initialise git fresh on your Mac (`.gitignore` is already set up here), create an **empty**
repo on GitHub (no README), then push:

```bash
cd ~/Desktop/Ovall26/newLES_cg
rm -rf .git                      # clear any partial repo state, start clean
git init && git add -A && git commit -m "Geometric SGD-plume LES (Oceananigans 0.109, CG solver)"
git branch -M main
git remote add origin git@github.com:<you>/iceplume_cg.git   # or the https URL
git push -u origin main
```

(Or with the GitHub CLI, after the `git init`/commit: `gh repo create iceplume_cg --private --source=. --push`.)

## 2. On Casper: clone and set up the environment (once)

```bash
# from a login node
cd /glade/work/$USER          # keep code + Julia depot on /glade/work (big, writable)
git clone git@github.com:<you>/iceplume_cg.git
cd iceplume_cg

# point JULIA at your Julia (install via juliaup, or use a build you have), then:
JULIA=/path/to/julia ./setup_casper.sh      # module loads + Pkg.instantiate + precompile
```

`setup_casper.sh` resolves the 0.109 environment and precompiles it. Do this on a login node
(or an interactive session) — it needs network for `Pkg`.

## 3. Submit the production runs

Edit `submit_pbs.sh` header if needed (`#PBS -A <account>`, walltime), then:

```bash
qsub -v CASE=vertical,JULIA=/path/to/julia submit_pbs.sh
qsub -v CASE=overcut,JULIA=/path/to/julia  submit_pbs.sh
```

Each job runs one case on one A100 with the GPU-scale domain (743 m × 4 km, 0.75 m fine to
375 m). It stops at 45 min model time or 11.5 h wall, and **auto-resumes from a checkpoint** if
re-submitted. Watch `logs/<case>.out`.

## 4. Get results back

Outputs + checkpoints land in **`/glade/work/$USER/LESplume_runs/`** (set by `OUTDIR` in
`submit_pbs.sh`, off the git repo). The 2-D slices and 15-min time-average are the light files;
the 3-D `*_fields.nc` are ~7 GB each (15-min cadence). Pull the light ones to your laptop and
plot (and copy anything you want to keep long-term to campaign `/glade/campaign/univ/uosc0035`):

```bash
# on your laptop, from the repo folder
scp 'derecho.hpc.ucar.edu:/glade/work/$USER/LESplume_runs/*_{midy,face,timeavg}.nc' output/
python3 plot_quicklook.py cg_overcut634 --dir output
python3 plot_quicklook.py cg_vertical   --dir output
```

## Notes / gotchas

- **Julia**: NCAR doesn't always provide a Julia module; installing via
  [juliaup](https://github.com/JuliaLang/juliaup) into `/glade/work/$USER` is simplest. Use
  Julia ≥ 1.10.
- **Depot on /glade/work**: `$HOME` has small quota; `setup_casper.sh` and `submit_pbs.sh`
  default `JULIA_DEPOT_PATH=/glade/work/$USER/.julia`.
- **Memory**: the default GPU domain is ~153 M cells → 80 GB A100. If you land on a 40 GB card
  or hit OOM, submit with a narrower fjord, e.g. add `--Ly=384` in `submit_pbs.sh` (~79 M).
- **First GPU compile** is slow (CUDA kernels); the `--pkgimages=no` flag in `submit_pbs.sh`
  avoids a known precompile issue on the cluster.
- **Reproducibility**: `Manifest.toml` is gitignored so laptop/cluster Julia versions don't
  clash. If you want an exact pin later, commit the cluster-resolved Manifest on a branch.
