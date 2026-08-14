#!/bin/bash
# Set up the Gaia DR4 completeness study environment on a cluster.
# Run this once on the login node, before submitting the array job.

set -e

# Octofitter v9 requires Julia >= 1.11.
module load julia/1.11.5
export JULIA_DEPOT_PATH="/scratch/$USER/julia_depot"
export DATADEPS_ALWAYS_ACCEPT=true

# Where the study runs. The scripts read and write their caches, results and
# plots relative to their own directory, so copy them here.
WORKDIR=/scratch/$USER/completeness_dr4
SRCDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$WORKDIR/logs"
cd "$WORKDIR"

cp "$SRCDIR"/*.jl .
cp "$SRCDIR"/submit.sh .

cat > Project.toml << 'PROJ'
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Octofitter = "daf3887e-d01a-44a1-9d7e-98f15c5d69c9"
Serialization = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
PROJ

# v9 is not registered yet; once it is, this is just `Pkg.add("Octofitter")`.
julia --project=. -e '
using Pkg
Pkg.add(url="https://github.com/sefffal/Octofitter.jl", rev="v2")
Pkg.instantiate()
Pkg.precompile()
'

# Cache the DR3 solution, the noise budget and the GOST scan forecast on the
# login node: compute nodes usually have no internet.
julia --project=. setup.jl

echo ""
echo "=== Setup complete ==="
echo "Smoke-test the whole workflow first (a few minutes):"
echo "  cd $WORKDIR && OCTO_COMPLETENESS_QUICK=1 julia --project=. run_local.jl"
echo "Then submit the array job:"
echo "  cd $WORKDIR && sbatch submit.sh"
