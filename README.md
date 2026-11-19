# Felodipine ER: PopPK QC

This branch holds only the PopPK part of the felodipine MBBE work, rewritten as six short steps
for QC. The PK-Sim model, the 2-D safe-space map, decks and documents stay on `main`.

**Question.** Do a PopPK model and the PBPK model (PK-Sim) reach the same bioequivalence (BE)
decisions for felodipine 10 mg ER, a highly variable BCS II drug?

**Approach.** PK-Sim simulated 27 ER formulations on a T50 × Weibull-shape grid; the reference is
T50 4.78 h, shape 0.764. Each formulation was run in the same 100 virtual subjects. We sample
those profiles like a BE study, fit the final PopPK model (M5) to each formulation, then run
identical BE code on the PK-Sim endpoints and on the PopPK endpoints.

## Files: read in this order

| # | File | What it does | Output (in `results/`) |
|---|---|---|---|
| 0 | `settings.jl` | Every number used: dose, sampling times, error, CVs, study sizes, seeds | — |
| 0 | `model_m5.jl` | The PopPK model (equations in the header) | — |
| 0 | `be_functions.jl` | Fast BE tests for the simulated trials (2x2 ABE, RSABE, ABEL), interpolation along an axis | — |
| 1 | `step1_build_dataset.jl` | PK-Sim Cmax/AUC per subject; samples the profiles at 18 times and adds 15% error | `step1_pbpk_endpoints.csv`, `step1_dataset.csv` |
| 2 | `step2_fit_poppk.jl` | Fits M5 to each formulation (FOCE), ~20–30 min per formulation | `fits/<T50>_<shape>/`, `step2_parameters.csv` |
| 3 | `step3_endpoints.jl` | Cmax and AUC 0–96 h per subject from the PopPK predictions (Pumas `NCA`); checks them against PK-Sim | `step3_poppk_endpoints.csv`, `step3_check_vs_pbpk.csv` |
| 4 | `step4_gmr.jl` | Test/reference GMR with 90% CI for all 27 formulations, PBPK vs PopPK (Pumas `pumas_be`) | `step4_gmr.csv`, `fig_step4_gmr.png` |
| 5 | `step5_power_safe_space.jl` | Power per formulation (ABE n = 54; RSABE/ABEL n = 24) and along the two axes | `step5_power_by_formulation.csv`, `step5_power_axes.csv`, figure |
| 6 | `step6_type1_error.jl` | BE boundary on each axis and the type I error there | `step6_boundaries.csv`, `step6_type1_error.csv`, `step6_tost_control.csv`, `step6_scaled_at_boundary.csv`, figures |
| 6b | `step6b_check_vs_pumas_be.jl` | Cross-check: our fast ABE, RSABE and ABEL tests vs `pumas_be` on the same simulated trials, at every boundary | `step6b_check_vs_pumas_be.csv` |

Inputs are in `data/pksim_exports/`: for each formulation, the PK-Sim population profiles
(`-Results.csv`) and PK-Sim's own NCA (`-PK-Analyses.csv`). The file name gives the formulation,
and a missing token means the reference value (for example, `Pop_b0.5` = T50 4.78 h, shape 0.5).

## The PopPK model (M5)

- **Release:** two Weibull components, fast and slow, mixed by the fraction FF and written with
  Pumas `@delay`, after a fixed 0.3113 h lag.
- **Absorption and disposition:** first-order absorption (ka) into a two-compartment model.
- **Variability:** on release, CL and Vc, with a proportional error.
- **Fitting:** one fit per formulation (100 subjects each).

M5 was chosen from six structures (M1–M6); the model testing log on `main` explains why.

## How to run

You need Julia 1.11 with Pumas. Step 2 uses `Pumas`, step 3 `NCA`, steps 4 and 6b `Bioequivalence`
(all part of the Pumas distribution). From this folder:

```bash
export JULIA="julia --project=."     # or your Pumas command, e.g. with a Pumas sysimage
./run_all.sh                         # all steps; step 2 takes ~10 h on one core
STEPS="step4_gmr" ./run_all.sh       # one step only
```

Step 2 can be split across cores by naming formulations, for example
`$JULIA step2_fit_poppk.jl 4.78_0.764 4.07_0.6`. Formulations already fitted are skipped. The
dense predictions (`fits/*/ipred.csv`, ~190 MB) are not committed; step 2 rewrites them.

## QC check: rerun vs the saved results

On 19 Nov 2026 all steps were rerun and compared with the saved M5 results
(`learning/Felodipine-BCSII/poppk/output_rich_m5/` and the `output_sweep_2d*` folders):

| Compared | Result |
|---|---|
| Dataset and PK-Sim endpoints (step 1) | identical |
| M5 fits, 27 formulations (step 2) | −2LL within 0.003 for all; parameters within 0.1% for 24. Three (T50 4.07 h at shape 0.9 and 1.0, 4.78 h at shape 1.2) differ by 1.5–4.7% at the same −2LL: a flat optimum |
| Cmax and AUC per subject (step 3) | within 0.04% |
| GMR and 90% CI (step 4) | within 0.05% |
| Power per formulation and along the axes (step 5) | within 0.1 percentage point (1 trial in 1,000) |
| Boundaries, type I error, TOST control, RSABE/ABEL at the boundary (step 6) | identical |

### Standard Pumas functions (19 Nov 2026)

- **Step 3: `NCA.cmax` and `NCA.auc`** (linear trapezoid, as before) replace our own code. On the
  same predictions they give identical Cmax and AUC (≤ 10⁻¹⁴).
- **Step 4: `pumas_be`** (2x2 crossover ANOVA) replaces our paired t-test. The reference row
  (reference vs itself) is dropped, since it has no CI.
- **Steps 5–6 keep the fast tests** in `be_functions.jl` for their hundreds of thousands of trials.
  **Step 6b** runs 2,000 trials at every BE boundary through both our test and `pumas_be`:

| Method | Same pass/fail decision | Largest pass-rate difference |
|---|---|---|
| ABE (2x2), Cmax and AUC | ≥ 99.45% of trials | 0.3 points |
| RSABE (FDA), `pumas_be` + `FDA_HighlyVariable` | ≥ 97.2% | 0.35 points |
| ABEL (EMA), EMA limits on `pumas_be`'s CI and CVwR | ≥ 97.65% | 1.85 points |

**Note on ABEL.** In this version of the Bioequivalence package, `pumas_be` with
`EMA_HighlyVariable` is not EMA ABEL. It runs the FDA RSABE method with θ = 0.76, without the EMA
widened limits or the 50% CV cap, and it prints a warning saying so. Step 6b therefore applies the
EMA rule (CVwR > 30%: CI within exp(±0.760·sWR), sWR capped at CVwR 50%; GMR within 0.80–1.25) to
`pumas_be`'s own estimates. The scaled methods agree slightly less than ABE because `pumas_be` fits
a mixed model to the replicate design, while our test uses within-subject differences.

## Key results (rich sampling, fasted)

PopPK/PBPK geometric-mean ratios over the 27 formulations: Cmax 0.976–1.019, AUC 0.989–1.004.

Formulations (of 27) with power ≥ 80%:

| Method | PBPK | PopPK | Same decision |
|---|---|---|---|
| ABE, Cmax + AUC (n = 54) | 6 | 5 | 26 / 27 |
| RSABE, Cmax (n = 24, replicate) | 12 | 10 | 23 / 27 |
| ABEL, Cmax (n = 24, replicate) | 9 | 8 | 26 / 27 |

Every disagreement goes the same way: PopPK fails a formulation that PBPK passes. PopPK is
slightly stricter, never more lenient.

BE boundaries (joint ABE; Cmax binds, AUC never leaves 0.80–1.25):

| Axis | Side | PBPK | PopPK | PBPK (true) pass rate at the PopPK boundary, n = 54 |
|---|---|---|---|---|
| T50 (h) | low | 3.25 | 3.27 | 5.6% |
| T50 (h) | high | 6.96 | 6.86 | 6.0% |
| shape | high | 1.46 | not reached | — |

A well-placed boundary gives 5%. On the T50 axis both PopPK boundaries sit within 0.1 h of PBPK's,
just inside it (true pass rate 5.6–6.0%, slightly strict). On the shape axis PopPK never reaches
the limit inside the simulated range (shape ≤ 1.5); PBPK reaches it at 1.46.

At the ABE Cmax boundary, scaled methods still pass often: RSABE 20–22%, ABEL 13–14% (n = 24), the
same for PBPK and PopPK. That is by design: they widen the limits for a highly variable drug.

## What to check

1. `step1`: dose 10 mg = 26,024 nmol; µmol/L → nM (×1000), µmol·min/L → nM·h (×1000/60);
   error applied as `true × (1 + 0.15·z)`.
2. `model_m5.jl`: the two-Weibull release, the lag, `bioav = 0` on Gut (the dose enters only
   through `@delay`), and units (h, nmol, L, nM).
3. `step2`: lag fixed; FOCE; starting values scale with the formulation's T50.
4. `step3`: PopPK/PBPK geometric-mean ratios close to 1 (`step3_check_vs_pbpk.csv`).
5. `be_functions.jl`: the trial simulation, RSABE (FDA) and ABEL (EMA) rules, and the
   interpolation along an axis; `step6b` shows the tests agree with `pumas_be`.
6. `step6`: the TOST control passes about 5% (reference × 0.80).

Not in this branch: the interpolated 2-D safe-space map, sparse sampling, the fed-state
sensitivity, and the PK-Sim model build.
