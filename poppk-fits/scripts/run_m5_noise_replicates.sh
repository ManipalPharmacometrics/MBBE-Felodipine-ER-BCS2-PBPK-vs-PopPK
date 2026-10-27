#!/bin/bash
REPO="$(cd "$(dirname "$0")/../.." && pwd)"   # repository root
# run_m5_noise_replicates.sh — is M5's shape-1.5 Cmax GMR shortfall (−3 % vs PBPK)
# noise or structure? Rebuild the dataset with two other residual-noise seeds and
# refit M5 on the arms that set the shape-axis BE limit (ref 0.764, 1.2, 1.5).
# Outputs output_m5_seed<seed>/; sentinel output/logs/m5_noise_replicates.DONE
J="julia -J${PUMAS_SO:-/opt/julia/lib/julia/pumas.so} --project="$REPO"/julia-env"
cd "$REPO"/poppk-fits
pids=()
for seed in 101 202; do
  (POPPK_OUT=$PWD/output_m5_seed$seed POPPK_SEED=$seed $J scripts/01_build_dataset.jl &&
   POPPK_OUT=$PWD/output_m5_seed$seed POPPK_REUSE=0 POPPK_ARMS=4.78_0.764,4.78_1.2,4.78_1.5 JULIA_NUM_THREADS=1 \
     $J scripts/02c_fit_poppk_m5.jl) > output/logs/m5_seed$seed.log 2>&1 &
  pids+=($!)
done
wait "${pids[@]}"
touch output/logs/m5_noise_replicates.DONE
