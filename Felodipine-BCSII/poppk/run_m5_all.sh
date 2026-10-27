#!/bin/bash
REPO="$(cd "$(dirname "$0")/../.." && pwd)"   # repository root
# run_m5_all.sh — final model M5 on all 27 arms (two shards, one per core), then
# collect to poppk/output_rich_m5, felo_02/03/05 (felo_04 BE compare run separately).
# Sentinel: poppk/run_m5_all.DONE (or .FAILED).
set -e
T="$REPO"/poppk-fits
F="$REPO"/Felodipine-BCSII
J="julia -J${PUMAS_SO:-/opt/julia/lib/julia/pumas.so} --project="$REPO"/julia-env"
trap 'touch $F/poppk/run_m5_all.FAILED' ERR
cd $T
echo "=== $(date -u +%H:%M) fits"
POPPK_SHARD=1/2 JULIA_NUM_THREADS=1 $J scripts/02c_fit_poppk_m5.jl > output/logs/02c_m5_shard1.log 2>&1 &
p1=$!
POPPK_SHARD=2/2 JULIA_NUM_THREADS=1 $J scripts/02c_fit_poppk_m5.jl > output/logs/02c_m5_shard2.log 2>&1 &
p2=$!
wait $p1 $p2
# rerun once unsharded to pick up any arm a shard missed (skips done arms)
$J scripts/02c_fit_poppk_m5.jl > output/logs/02c_m5_sweep.log 2>&1
echo "=== $(date -u +%H:%M) collect: $(ls output_m5/arms/*/nca.csv | wc -l) arms"
POPPK_OUT=$T/output_m5 LEARN_DST=$F/poppk/output_rich_m5 $J scripts/03_collect.jl
cd "$REPO"/programs
run() { echo "=== $(date -u +%H:%M) $*"; julia --project=. "$@"; }
run $F/scripts/felo_02_be_sweep.jl $F/poppk/output_rich_m5/pk_analyses_format $F/poppk/output_sweep_2d_poppk_m5
run $F/scripts/felo_03_safespace_2d.jl $F/poppk/output_sweep_2d_poppk_m5
run $F/scripts/felo_03_safespace_2d.jl $F/output_sweep_2d $F/poppk/output_sweep_2d_poppk_m5 PopPK-M5
run $F/scripts/felo_05_poppk_diagnostics.jl $F/poppk/output_rich_m5
echo "=== $(date -u +%H:%M) ALL DONE"
touch $F/poppk/run_m5_all.DONE
