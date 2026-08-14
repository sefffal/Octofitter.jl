#!/bin/bash
#SBATCH --job-name=dr4_compl
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-720
#SBATCH --output=logs/trial_%a_%j.out
#SBATCH --error=logs/trial_%a_%j.err
#SBATCH --account=def-cmarois-ab

# Gaia DR4 completeness mapping — SLURM array job.
#
# The array range must match the grid in common.jl:
#   12 masses × 12 separations × 5 trials = 720 jobs.
# `setup.jl` prints the number to use.

# Octofitter v9 requires Julia >= 1.11.
module load julia/1.11.5
export JULIA_DEPOT_PATH="/scratch/$USER/julia_depot"
export JULIA_NUM_THREADS=4
export DATADEPS_ALWAYS_ACCEPT=true

cd /scratch/$USER/completeness_dr4

julia --project=. --threads=4 completeness_trial.jl

echo "Trial $SLURM_ARRAY_TASK_ID completed at $(date)"
