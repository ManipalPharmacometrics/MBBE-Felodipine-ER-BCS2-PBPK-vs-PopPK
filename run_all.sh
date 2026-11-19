#!/usr/bin/env bash
# Run steps 1–6 in order. Each step's console output goes to results/logs/.
# Set JULIA to your Pumas-enabled Julia command first, for example:
#   export JULIA="julia --project=."
set -euo pipefail
cd "$(dirname "$0")"
JULIA=${JULIA:-"julia --project=."}
mkdir -p results/logs
for step in ${STEPS:-step1_build_dataset step2_fit_poppk step3_endpoints step4_gmr step5_power_safe_space step6_type1_error step6b_check_vs_pumas_be}; do
    echo "$(date +%H:%M:%S)  $step"
    $JULIA "$step.jl" > "results/logs/$step.log" 2>&1
done
echo "$(date +%H:%M:%S)  done"
