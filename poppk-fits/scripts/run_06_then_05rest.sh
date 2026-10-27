#!/bin/bash
REPO="$(cd "$(dirname "$0")/../.." && pwd)"   # repository root
# 2026-10-24: structural variants first, then the rest of the no-lag @delay fits.
cd "$REPO"/poppk-fits
J="julia -J${PUMAS_SO:-/opt/julia/lib/julia/pumas.so} --project="$REPO"/julia-env"
$J scripts/06_structure_fits.jl > output/logs/06_structure_fits.log 2>&1
MC_MODELS=M3d0,M2d0 $J scripts/05_delay_model_fits.jl > output/logs/05_delay_model_fits_rest.log 2>&1
touch output/logs/run_06_then_05rest.DONE
