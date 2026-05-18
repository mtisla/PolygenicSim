#!/bin/bash
#SBATCH -J recap_ism_s3
#SBATCH --account=rogersa
#SBATCH --partition=kingspeak
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8gb
#SBATCH --output=out/s3/s3.slurm.out
#SBATCH -e out/s3/s3.slurm.err
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
SEED=3
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
