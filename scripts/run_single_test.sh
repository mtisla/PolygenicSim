#!/usr/bin/env bash
# Run exactly one @testset from test/runtests.jl in isolation.
#
# Every @testset in test/runtests.jl is self-contained: the only state it
# shares with the rest of the file is the `using ...` block and
# `const PS = PolygenicSim` at the top (verified — no shared helper
# functions or global mutable state across testsets as of v0.21.0). That
# means a single testset can be extracted by name and run standalone with
# zero edits to the source file: this script greps out the block between
# `@testset "<name>" begin` and its closing `end` and re-runs it under a
# fresh `using` preamble.
#
# Usage:
#   scripts/run_single_test.sh --list                       # print every testset name + line
#   scripts/run_single_test.sh "Test 4 — Haldane recomb fraction"
#   JULIA_NUM_THREADS=8 scripts/run_single_test.sh "Test 9 — dense ≡ packed (bit-identical)"
#
# JULIA_NUM_THREADS defaults to 4 to match the project convention (the sim
# is deterministic per thread count, not across thread counts — see
# README "Reproducibility"). Override by exporting it before calling this
# script if you need a different pinned count.

set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

TESTFILE="test/runtests.jl"

if [[ $# -eq 0 || "$1" == "--list" || "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Available testsets in ${TESTFILE} (name — line):" >&2
  grep -n '@testset "' "$TESTFILE" | tail -n +2 | \
    sed -E 's/^([0-9]+):@testset "(.*)" begin/  L\1\t\2/'
  echo "" >&2
  echo "Usage: $0 \"<exact testset name>\"" >&2
  exit 0
fi

NAME="$1"
OPENER="@testset \"${NAME}\" begin"

BLOCK=$(awk -v t="$OPENER" '
  index($0, t) == 1 { p = 1 }
  p { print }
  p && $0 == "end" { exit }
' "$TESTFILE")

if [[ -z "$BLOCK" ]]; then
  echo "error: no testset named '${NAME}' found in ${TESTFILE}" >&2
  echo "  run '$0 --list' to see valid names (must match exactly, including em-dashes)" >&2
  exit 1
fi

TMP=$(mktemp /tmp/polygenicsim_single_test.XXXXXX.jl)
trap 'rm -f "$TMP"' EXIT
printf '%s\n' "$BLOCK" > "$TMP"

export JULIA_NUM_THREADS="${JULIA_NUM_THREADS:-4}"
echo ">> running '${NAME}' at JULIA_NUM_THREADS=${JULIA_NUM_THREADS}" >&2

julia --project=. -e '
using Test, Random, Statistics, StatsBase, Distributions, TOML, PolygenicSim
const PS = PolygenicSim
include(ARGS[1])
' "$TMP"
