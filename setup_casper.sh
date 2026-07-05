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

# --- modules (match submit_pbs.sh) ---
module --force purge || true
module load ncarenv/23.10 || true
module load cuda || true

# --- Julia: point JULIA at your build/module; keep the depot on /glade/work (big, writable) ---
JULIA="${JULIA:-julia}"                                   # e.g. /glade/u/home/$USER/bin/julia-1.10.x/bin/julia
export JULIA_DEPOT_PATH="${JULIA_DEPOT_PATH:-/glade/work/$USER/.julia}"
mkdir -p "$JULIA_DEPOT_PATH"
echo "Julia:  $($JULIA --version)"
echo "Depot:  $JULIA_DEPOT_PATH"

# --- resolve + build the env (no Manifest is shipped, so this resolves 0.109.x fresh) ---
$JULIA --project -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

echo
echo "Environment ready. Submit the production runs (edit walltime/account in submit_pbs.sh first):"
echo "    qsub -v CASE=vertical submit_pbs.sh"
echo "    qsub -v CASE=overcut  submit_pbs.sh"
echo "Tip: a quick GPU sanity check before the big job —"
echo "    $JULIA --project iceplume.jl --arch=gpu --face_angle=63.4 --Ly=384 --Lx=1500 --stop_time=2 --simname=gpucheck"
