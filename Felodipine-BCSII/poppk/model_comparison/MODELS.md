# Felodipine ER PopPK: candidate models

Status as of 2026-10-25. Pumas sources live in `poppk-fits/src/`, with copies in `../src/`.

**Final model: M5** (dual Weibull in the `@delay` form, fixed lag → ka → 2-cmt, per arm), fitted on all 27 arms. **Fallback: M3.** Equations and code for every model: `model_equations/felodipine_poppk_models.pdf` (or `.html`).

## What every model shares

**Data.** Each arm has 100 PK-Sim subjects with 18 samples per subject (0.5–96 h) and 15% proportional noise. The dose is 10 mg = 26 024 nmol. Concentrations are in nmol/L, time in h, volumes in L, and clearances in L/h.

**Release lag.** t_lag = 0.3113 h, held fixed at the PK-Sim product value.

**Weibull release.** Let τ = t − t_lag. The hazard is

    h(τ) = (γ/α) · (τ/α)^(γ−1)   for τ > 0, else 0
    α = T50 / (ln 2)^(1/γ)        (M3/M4 estimate T50; M1/M2 estimate the scale α)

The release clock is floored at τ ≥ 1e-8 h in M1/M2 and τ ≥ 0.01 h in M3/M4. Without the larger floor, γ < 1 gives an Inf likelihood at the start values.

**Residual error.** y = C · (1 + ε), ε ~ N(0, σ²), so the residual SD is σ·C + 1e-6.

**IIV.** Log-normal with diagonal Ω: P_i = θ_P · exp(η_P), η ~ N(0, Ω).

**Estimation.** FOCE from two starts (listed per model); the fit with the higher log-likelihood is kept.

All volumes and clearances are apparent (/F).

---

## M1: 1-cmt + Weibull release (v1 baseline, 2026-10-16)

`poppk_model_felo_1cmt.jl`

    dDepot/dt   = −h(τ) · Depot
    dCentral/dt =  h(τ) · Depot − (CL/Vc) · Central
    C = Central / Vc

- **Fixed effects:** α, γ, CL, Vc (4).
- **IIV:** α, CL, Vc.
- **Parameter count:** 8 (θ 4, ω 3, σ 1).
- **Starts:**
  - release-limited: α 13, γ 1.0, Vc 250.
  - disposition-limited: α 1, γ 1.3, Vc 4000.
- **Result:**
  - Two basins, and which one wins depends on the arm.
  - In 23 of 27 arms the fit has a fast "release" (implied T50 ≈ 0.3–2.4 h) with Vc ≈ 4000 L.
  - Cmax/PBPK runs 0.87–1.135 across arms, which biases the T/R GMR by up to 20%.

## M2: 2-cmt + Weibull release (tried 2026-10-16, rejected)

`poppk_model_felo.jl`

    dDepot/dt   = −h(τ) · Depot
    dCentral/dt =  h(τ) · Depot − (CL/Vc) · Central − (Q/Vc) · Central + (Q/Vp) · Periph
    dPeriph/dt  =  (Q/Vc) · Central − (Q/Vp) · Periph
    C = Central / Vc

- **Fixed effects:** α, γ, CL, Vc, Q, Vp (6).
- **IIV:** α, CL, Vc.
- **Parameter count:** 10.
- **Starts:**
  - release-limited: α 7, γ 0.9, Vc 250, Q 200, Vp 1500.
  - disposition-limited: α 1, γ 1.3, Vc 2000, Q 300, Vp 3000.
- **Result on 10-16:** Q or Vp ran to their bounds from either start, so the peripheral compartment is not identified. Refitted here on the diagnostic arms so it can sit in the table.

## M3: Weibull release → first-order absorption → 2-cmt, per arm (v2, 2026-10-23)

`poppk_model_felo_v2.jl` → `poppk_felo_wka2`

    dDepot/dt   = −h(τ) · Depot
    dGut/dt     =  h(τ) · Depot − ka · Gut
    dCentral/dt =  ka · Gut − ((CL + Q)/Vc) · Central + (Q/Vp) · Periph
    dPeriph/dt  =  (Q/Vc) · Central − (Q/Vp) · Periph
    C = Central / Vc

- **Fixed effects:** T50, γ, ka, CL, Vc, Q, Vp (7).
- **IIV:** T50, CL, Vc.
- **Parameter count:** 11.
- **Starts** (release starts at the arm's PK-Sim T50 and shape):
  - release-limited: ka 2, Vc 200, Q 200, Vp 1500.
  - absorption-limited: ka 0.3, Vc 500, Q 100, Vp 1000.
- **Rationale:** the gut depot separates release from appearance in plasma. That gives the PBPK profile's fast rise followed by a plateau, without forcing release to be fast.

## M4: M3 paired, reference + test in the same subjects (prototype, 2026-10-23)

`poppk_model_felo_v2.jl` → `poppk_felo_pair`

The reference is dosed at t = 0 into Depot1. The test is dosed at t = 200 h into Depot2, with τ₂ = t − 200 − t_lag. The two depots release with their own Weibull parameters.

    dDepot1/dt  = −h_R(τ)  · Depot1         h_R: T50_R, γ_R
    dDepot2/dt  = −h_T(τ₂) · Depot2         h_T: T50_T, γ_T
    dGut/dt     =  h_R · Depot1 + h_T · Depot2 − ka · Gut
    dCentral/dt =  ka · Gut − ((CL + Q)/Vc) · Central + (Q/Vp) · Periph
    dPeriph/dt  =  (Q/Vc) · Central − (Q/Vp) · Periph

- **Shared between periods:** ka, CL, Vc, Q, Vp, and each subject's η.
- **IIV:** η_T50 (applied to both T50_R and T50_T), CL, Vc.
- **Fixed effects:** T50_R, γ_R, T50_T, γ_T, ka, CL, Vc, Q, Vp (9).
- **Parameter count:** 13.
- **Starts:** as M3.
- **AIC:** M4 is fitted to two arms of data, so its AIC can't be compared with the per-arm models. Compare its LL with the sum of the M3 LLs for the reference arm and the test arm.
- **Prototype result:** LL is far worse than the two M3 fits combined, σ rises to 0.18–0.25, and the 1.5 Cmax GMR is 1.15 (PBPK 1.27).
  - Why: in the PBPK model, release rate changes where along the gut the drug is absorbed and how much CYP3A4 gut-wall extraction it sees. The disposition seen from the gut is therefore not shared across formulations.

## M3d: M3 in the `@delay` form (2026-10-24)

`poppk_models_felo_delay.jl` → `poppk_m3d` (also M1d, M2d)

The same equations as M3, but the release goes in as `@delay(Tlag + truncated(Weibull(γ, α), 1e-4, Inf), Gut)`, with the dose record into Gut at bioavailability 0. This is the `@delay` form.

- `@dosecontrol` lags are ignored by `@delay`, so the lag sits inside the delay distribution.
- Without the truncation, shape < 1 gives an Inf likelihood. The truncation removes ≤ 0.4% of the release.
- **Result:** reproduces M3 on the 5 diagnostic arms (best LL within 3, identical Tmax and curves). Without the lag, the fits and curves are much worse. Estimating the lag fails in this form; in the ODE form (M3Lo) it lands at about 0.3 h and the fits get worse, so the lag stays fixed.

## M5: dual Weibull → first-order absorption → 2-cmt (FINAL, 2026-10-24/25)

`poppk_models_felo_structure.jl` → `poppk_m5`

    I(t)        = Dose · [FF · f_W(τ; γa, αa) + (1 − FF) · f_W(τ; γb, αb)]   (`@delay`, truncated at 1e-4 h)
    dGut/dt     = I(t) − ka · Gut
    dCentral/dt = ka · Gut − ((CL + Q)/Vc) · Central + (Q/Vp) · Periph
    dPeriph/dt  = (Q/Vc) · Central − (Q/Vp) · Periph

- **Fixed effects:** T50a, γa, T50b, γb, FF, ka, CL, Vc, Q, Vp (10).
- **IIV:** one η on both T50s, plus CL and Vc.
- **Parameter count:** 14.
- **Starts:** release-limited (T50a 0.5T, T50b 3T, FF 0.6, ka 2) and absorption-limited (T50a 0.3T, T50b 2T, γ 1, FF 0.3, ka 0.5, Vc 500). The 27-arm run (`02c_fit_poppk_m5.jl`, `output_m5/`) used the absorption-limited start only; on the 7 arms fitted from both starts it was never worse.
- **Result, all 27 arms vs M3:**
  - AIC better on 27/27 (median −48).
  - Cmax/PBPK 0.976–1.019 (M3 0.976–1.043).
  - GMR error max 3.3% (M3 3.5%).
  - Tmax bias median of arms 0.20 h (M3 0.45 h).
  - Safe-space limits no closer to PBPK: mean |Δ| 0.09 h, shifted −0.09 h (M3 0.08 h, +0.03 h), because the reference-arm Cmax is +0.9%.
- **Caveats:**
  - The two components are descriptive (FF 0.02–0.97 across arms).
  - The mixture's T50 sits at ≈ 8–9 h vs PK-Sim 4–5 h, as for M3.
  - FF and Q/Vp are poorly determined at the reference arm.
  - One arm, 6.31/1.2, still peaks 0.85 h early.

## M6: Weibull → 3 transit compartments → 2-cmt (tested, not selected)

`poppk_models_felo_structure.jl` → `poppk_m6`. The transit chain (one ktr) replaces ka; 11 parameters. It fixes the shape-1.5 curve but is no better than M3 on the other diagnostic arms, and 2 of its 10 fits failed.

---

## Files in this folder

| File | Contents |
|---|---|
| `all27_M1_vs_M3.csv` | Production fits for all 27 arms, M1 vs M3: LL, AIC, BIC, σ, Cmax/AUC agreement with PBPK, and model vs PBPK T/R GMR |
| `diagnostic_fits.csv` | M1–M4 on 5 diagnostic arms (ref, 4.78/0.5, 4.78/1.5, 2.39/0.764, 9.56/0.764). One row per start: status, LL, AIC, BIC, σ, agreement, GMR error, ΔAIC within arm |
| `diagnostic_params.csv` | Every parameter estimate, one row per model × arm × start × parameter |
| `diagnostic_best.csv` | Best start per model × arm, with parameters as columns |
| `prototype_runs_2026-10-23.csv` | The first M3/M4 prototype fits, recovered from their logs |
| `fig_all27_M1_vs_M3.png` | Cmax agreement and GMR error by arm |
| `delay_vs_ode_fits.csv`, `delay_vs_ode_params.csv` | M1–M3 hazard ODE vs the `@delay` form, with lag (d) and without (d0) |
| `curve_tmax/` | Tmax and mean curve vs noise-free PK-Sim on the diagnostic arms, every model (`felo_08`) |
| `m5_vs_m3_all27/` | Final M5 vs M3 on all 27 arms: per-arm LL/AIC, agreement, GMR, Tmax; safe-space limits vs PBPK (`felo_09`) |
| `model_equations/` | Equations + Pumas code for M1–M6 and the model-selection summary (PDF/HTML; `build_model_equations.py`) |

The per-fit raw outputs for the diagnostic arms are in `poppk-fits/output_modelcomp/<model>/<arm>/<start>/`:

- `coef.csv`, `se.csv`, `shrinkage.csv`;
- `endpoints.csv` (per subject);
- `ipred_thin.csv`;
- `fit.jls` (serialized fit; reload in the julia-env environment).

To rebuild the tables, run from `programs`:

    julia --project=. ../Felodipine-BCSII/scripts/felo_07_model_comparison.jl
    julia --project=. ../Felodipine-BCSII/scripts/felo_08_curve_tmax_check.jl
    julia --project=. ../Felodipine-BCSII/scripts/felo_09_m5_vs_m3_all27.jl

The M5 27-arm results are in `../output_rich_m5/` (sweep `../output_sweep_2d_poppk_m5/`). The M3 ones are in `../output_rich_v2/`.
