#!/bin/bash -l
# setup_casper.sh — one-time environment setup on Casper (NCAR), run from the repo directory.
# Resolves + precompiles the Oceananigans 0.109 Julia environment for the GPU runs.
#
#   git clone <your-repo-url> iceplume_cg && cd iceplume_cg
#   JULIA=/path/to/julia ./setup_casper.sh
#   qsub -v CASE=vertical submit_pbs.sh
#   qsub -v CASE=overcut  submit_pbs.sh
set -e
cd "$(dirname "$0")"

# NOTE: instantiating Julia packages needs NO HPC modules (CUDA.jl bundles its own toolkit; the
# GPU driver is only needed at run time), so this script just needs a working `julia`.
# If `which julia` is empty, get one first:
#   module spider julia                         # is there an NCAR Julia module? if so, module load it
#   # else install juliaup (~1-2 GB, fine on the 100 GB home quota):
#   curl -fsSL https://install.julialang.org | sh -s -- --yes && source ~/.bashrc

# --- Julia (packages go to the default depot ~/.julia unless you export JULIA_DEPOT_PATH) ---
JULIA="${JULIA:-julia}"                                   # or a full path to your julia binary
command -v "$JULIA" >/dev/null 2>&1 || { echo "ERROR: '$JULIA' not found — install Julia (see notes above) then re-run."; exit 1; }
echo "Julia:  $($JULIA --version)   depot: ${JULIA_DEPOT_PATH:-$HOME/.julia}"

# --- resolve + build the env (no Manifest is shipped, so this resolves 0.109.x fresh) ---
$JULIA --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo
echo "Environment ready. Submit the production runs (edit walltime/account in submit_pbs.sh first):"
echo "    qsub -v CASE=vertical submit_pbs.sh"
echo "    qsub -v CASE=overcut  submit_pbs.sh"
echo "Tip: a quick GPU sanity check before the big job —"
echo "    $JULIA --project iceplume.jl --arch=gpu --face_angle=63.4 --Ly=384 --Lx=1500 --stop_time=2 --simname=gpucheck"
