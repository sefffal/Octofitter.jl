#!/bin/bash
# Setup the completeness study environment on fir.
# Run this once on the login node before submitting jobs.

set -e

module load julia/1.10.10
export JULIA_DEPOT_PATH="/scratch/wthompso/julia_depot"
export DATADEPS_ALWAYS_ACCEPT=true

WORKDIR=/scratch/wthompso/completeness_dr4

mkdir -p "$WORKDIR/results"
mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"

# Copy scripts from repo
cp ~/octo-3/Octofitter-v8/examples/completeness_dr4/*.jl .
cp ~/octo-3/Octofitter-v8/examples/completeness_dr4/submit.sh .

# Create Project.toml
cat > Project.toml << 'PROJ'
[deps]
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Octofitter = "daf3887e-d01a-44a1-9d7e-98f15c5d69c9"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
PROJ

# Add Octofitter from the completeness-map branch
julia --project=. -e '
using Pkg
Pkg.add(url="https://github.com/sefffal/Octofitter.jl", rev="completeness-map")
Pkg.add("CairoMakie")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.instantiate()
Pkg.precompile()
'

# Pre-cache GOST scan law on login node (compute nodes may not have internet)
julia --project=. setup.jl

echo ""
echo "=== Setup complete ==="
echo "Submit the array job with:"
echo "  cd $WORKDIR && sbatch submit.sh"
