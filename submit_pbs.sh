#!/bin/bash -l
#PBS -A UGIT0046
#PBS -N iceplume_cg
#PBS -o logs/iceplume_cg.log
#PBS -j oe
#PBS -l walltime=11:59:00
#PBS -q casper
#PBS -l select=1:ncpus=4:mem=60gb:ngpus=1:gpu_type=a100
#PBS -M kenzhao@unc.edu
# mem=60gb is HOST RAM (not GPU): the GPU→host staging for 3-D NetCDF output and, especially,
# the full-state checkpoint of a 153 M-cell run needs ~15-30 GB. The default (~10 gb) would OOM
# at the first checkpoint. The A100-80GB nodes have plenty of host RAM; raise to 120gb if needed.
#PBS -m ae
#PBS -r n

# PBS starts the job in $HOME, so cd to the submission directory (where you ran qsub) — that's
# where iceplume.jl, Project.toml, and logs/ live. logs/ must exist because PBS also writes its
# -o/-e files there; a missing logs/ is what makes the job requeue (R->Q).
cd "$PBS_O_WORKDIR" || exit 1
mkdir -p logs

# Geometric-terminus SGD-plume LES (Oceananigans 0.109, CG solver). One job = one case:
#     qsub -v CASE=vertical submit_pbs.sh
#     qsub -v CASE=overcut  submit_pbs.sh     (default)
#
# GPU-ONLY domain (the CPU-test defaults in iceplume.jl are left small on purpose):
#   743 m wide × 4 km long, 0.75 m in the 375 m nearest the terminus, then the along-fjord
#   grid stretches to 18.6 m out to 4 km. ~153 M cells — plan on an 80 GB A100 (a 40 GB card
#   will be tight). Each 3-D snapshot is ~7 GB, so 3-D output is written sparingly (15 min);
#   the 2-D face/mid-y slices and the 15-min time-average carry the detail.
#   If it's too heavy: drop to --Ly=384 (~79 M cells).

# No HPC modules are required: CUDA.jl bundles its own CUDA toolkit and uses the GPU node's
# driver (libcuda). If CUDA.jl ever can't find the driver, `module load cuda` (first check the
# exact name with `module avail cuda`) — but do NOT force-purge / load a fixed ncarenv version,
# which on current Casper pulls a broken openmpi.

# NOTE: this is a DIFFERENT (0.109) environment from ../newLES. Run ./setup_casper.sh once on
# the cluster first (instantiates the env; no Manifest is shipped — it resolves fresh).
# Point JULIA at your Julia binary (juliaup: $HOME/.juliaup/bin/julia). The package depot
# defaults to ~/.julia — the SAME one setup_casper.sh instantiated, so packages are found.
# Override via `qsub -v CASE=overcut,JULIA=/path/to/julia submit_pbs.sh`.
JULIA="${JULIA:-julia}"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES ; Julia: $($JULIA --version 2>/dev/null)"

CASE=${CASE:-overcut}
DOMAIN="--Ly=743 --Lx=4000 --fine_x=375 --dz=0.75 --dx_max=18.6"     # GPU-only domain
# Outputs + checkpoints go here (a new folder on /glade/work, OFF the git repo). Override with
# `qsub -v CASE=overcut,OUTDIR=/glade/derecho/scratch/$USER/LESplume submit_pbs.sh` for scratch.
OUTDIR="${OUTDIR:-/glade/work/$USER/LESplume_runs}"
mkdir -p "$OUTDIR"
if [ "$CASE" = "vertical" ]; then
    ARGS="--simname=cg_vertical --face_angle=90"
else
    ARGS="--simname=cg_overcut634 --terminus=overcut --face_angle=63.4"   # 2:1 headline overcut
fi

# Discharge 150 m3/s, 24x6 m point-plume outlet, grounding line 150 m (2024 conditions);
# 45 min run, 15-min time-average. Auto-resumes from a checkpoint if the job is re-submitted.
time $JULIA --project --pkgimages=no iceplume.jl \
    $ARGS $DOMAIN --discharge=150 --outlet_w=24 --outlet_h=6 --Lz=150 \
    --stop_time=45 --output_interval=900 --wall_time_limit=11.5 --outdir="$OUTDIR" \
    2>&1 | tee logs/${CASE}.out

qstat -f $PBS_JOBID >> logs/iceplume_cg.log
