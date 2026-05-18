#!/bin/bash
# Run AFTER all 3 seeds have finished. Produces the median-of-3 table.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
OUT_DIR="$SCRIPT_DIR/out"

# Sanity check.
for s in 1 2 3; do
    for f in h2.txt oracle.init.tsv oracle.settled.tsv oracle.1.0_thalf.tsv oracle.2.0_thalf.tsv; do
        if [ ! -f "$OUT_DIR/s$s/s$s.$f" ]; then
            echo "MISSING: out/s$s/s$s.$f"
            echo "Run seed $s before aggregating."
            exit 1
        fi
    done
done

julia --project="$REPO_ROOT" "$SCRIPT_DIR/aggregate.jl" "$OUT_DIR"
echo "Aggregated table at: $OUT_DIR/aggregated.tables.txt"
