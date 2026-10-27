#!/bin/bash
REPO="$(cd "$(dirname "$0")/../.." && pwd)"   # repository root
# run_analysis_v2.sh — BE sweep, 2-D safe space, PBPK-vs-PopPK compare + Type I,
# diagnostics and the model-comparison tables on the v2 PopPK (M3) results.
# Waits for the Pumas model-comparison fits (04) so the 2 cores are not shared.
set -e
while pgrep -f 04_model_comparison_fits.jl >/dev/null; do sleep 60; done
cd "$REPO"/programs
F=../Felodipine-BCSII
run() { echo "=== $(date -u +%H:%M) $*"; julia --project=. "$@"; }
run $F/scripts/felo_02_be_sweep.jl $F/poppk/output_rich_v2/pk_analyses_format $F/poppk/output_sweep_2d_poppk_v2
run $F/scripts/felo_03_safespace_2d.jl $F/poppk/output_sweep_2d_poppk_v2
run $F/scripts/felo_03_safespace_2d.jl $F/output_sweep_2d $F/poppk/output_sweep_2d_poppk_v2 PopPK-v2
run $F/scripts/felo_04_pbpk_poppk_compare.jl $F/poppk/output_rich_v2
run $F/scripts/felo_05_poppk_diagnostics.jl $F/poppk/output_rich_v2
run $F/scripts/felo_07_model_comparison.jl
echo "=== $(date -u +%H:%M) ALL DONE"
