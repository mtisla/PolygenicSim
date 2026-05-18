#!/bin/bash
#SBATCH -J recap_ism_s1
#SBATCH --account=rogersa
#SBATCH --partition=kingspeak
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --output=s1.slurm-%j.out
#SBATCH -e s1.slurm-%j.err
set -euo pipefail

# Load Julia on the compute node (sbatch jobs don't inherit login modules).
module load julia/1.11.1

# Resolve paths via SLURM_SUBMIT_DIR (the directory where `sbatch` was run).
# BASH_SOURCE points to SLURM's spool copy of the script on kingspeak,
# which is not writable by the user.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SEED=1
OUT_DIR="$SCRIPT_DIR/out/s${SEED}"
OUT_PREFIX="$OUT_DIR/s${SEED}"
mkdir -p "$OUT_DIR"

echo "[$(date)] seed=$SEED start"
export JULIA_NUM_THREADS=2
julia --project="$REPO_ROOT" "$SCRIPT_DIR/sim_seed.jl" "$SEED" "$OUT_PREFIX"
echo "[$(date)] sim done; formatting tables…"
julia --project="$REPO_ROOT" "$SCRIPT_DIR/format_seed.jl" "$OUT_PREFIX"
echo "[$(date)] seed=$SEED finished"
echo "Pretty table at: $OUT_PREFIX.tables.txt"
