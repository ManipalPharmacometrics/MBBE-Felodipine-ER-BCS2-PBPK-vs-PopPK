# build_model_equations.py — write felodipine_poppk_models.qmd (equations + the
# Pumas model code pasted verbatim from the fitted sources), 2026-10-23.
# 2026-10-25: added the `@delay` form (M3d), the structural variants M5/M6 and the
# final-model choice (M5, M3 fallback) with the 27-arm results.
# Render: quarto render felodipine_poppk_models.qmd   (-> .html and .pdf)
import pathlib

HERE = pathlib.Path(__file__).resolve().parent
SRC = (pathlib.Path(__file__).resolve().parent / "../../../../poppk-fits/src").resolve()
SCRIPTS = SRC.parent / "scripts"


def code(fname, start=None, end=None):
    """Return the source between the line containing `start` and the line containing `end` (inclusive)."""
    lines = (SRC / fname).read_text().splitlines()
    i = 0 if start is None else next(k for k, l in enumerate(lines) if start in l)
    j = len(lines) - 1 if end is None else next(k for k, l in enumerate(lines) if k >= i and l.startswith(end))
    return "\n".join(lines[i:j + 1])


M1 = code("poppk_model_felo_1cmt.jl", "poppk_model_felo = @model", "    \"disposition_limited\"")
M2 = code("poppk_model_felo.jl", "poppk_model_felo = @model", "end")
M3 = code("poppk_model_felo_v2.jl", "poppk_felo_wka2 = @model", "end")
M4 = code("poppk_model_felo_v2.jl", "poppk_felo_pair = @model", "end")
INITS = code("poppk_model_felo_v2.jl", "const INITS_V2", "    \"absorption_limited\"")
TP2 = code("poppk_model_felo_v2.jl", "const TP2", "const TP2")
M3D = code("poppk_models_felo_delay.jl", "poppk_m3d = @model", "end")
M5 = code("poppk_models_felo_structure.jl", "poppk_m5 = @model", "end")
M6 = code("poppk_models_felo_structure.jl", "poppk_m6 = @model", "end")


def script_block(fname, start, stop):
    """Lines of a script from the one containing `start` to the one containing `stop`."""
    lines = (SCRIPTS / fname).read_text().splitlines()
    i = next(k for k, l in enumerate(lines) if start in l)
    j = next(k for k, l in enumerate(lines) if k >= i and stop in l)
    return "\n".join(lines[i:j + 1])


M5STARTS = script_block("02c_fit_poppk_m5.jl", "starts(t, s) = [", "tvvp = 1000.0, tvlag = 0.3113)]")

doc = rf"""---
title: "Felodipine ER PopPK: candidate model structures"
subtitle: "M1–M6 equations and Pumas model code; final model M5 (fallback M3)"
date: 2026-10-25
format:
  html:
    toc: true
    embed-resources: true
    html-math-method: mathjax
    code-overflow: wrap
  typst:
    toc: true
    papersize: a4
    margin:
      x: 2cm
      y: 2cm
---

Source files: `poppk-fits/src/` (copies in
`Felodipine-BCSII/poppk/src/`). The code blocks below are pasted verbatim
from those files by `build_model_equations.py`.

**Final model (2026-10-25): M5**, dual Weibull release in the `@delay` form with the
fixed release lag → first-order absorption → 2-compartment, fitted per arm on all 27
arms. **M3** stays as the fallback. The results that led to this choice are summarised in
[Model selection](#model-selection) at the end.

# Shared elements

**Data.** 100 PK-Sim subjects per arm, 18 samples each (0.5–96 h), 15 % proportional noise.
Dose 10 mg = 26 024 nmol; units h, nmol, L, nmol/L.

**Release lag.** $t_{{\text{{lag}}}} = 0.3113$ h (fixed, PK-Sim product value); $\tau = t - t_{{\text{{lag}}}}$.

**Weibull input (`@delay` format).** Using the Pumas built-in delay input (`@delay(Weibull(SS, TD), Central)`), the dose enters
as $\text{{Dose}}$ times a Weibull density with shape $SS$ and scale $TD$ (Distributions.jl order
`Weibull(shape, scale)`), started at the lag:

$$
I(t) = \text{{Dose}}\cdot f_W(t - t_{{\text{{lag}}}};\, SS, TD), \qquad
f_W(\tau;\, SS, TD) =
\begin{{cases}}
  \dfrac{{SS}}{{TD}}\left(\dfrac{{\tau}}{{TD}}\right)^{{SS-1}} \exp\left[-\left(\dfrac{{\tau}}{{TD}}\right)^{{SS}}\right], & \tau > 0\\[6pt]
  0, & \tau \le 0
\end{{cases}}
$$

Cumulative fraction released/absorbed (`ABSO`):

$$
\text{{ABSO}}(t) = 1 - \exp\left[-\left(\frac{{t - t_{{\text{{lag}}}}}}{{TD}}\right)^{{SS}}\right],
\qquad
T_{{50}} = TD\,(\ln 2)^{{1/SS}} \;\Longleftrightarrow\; TD = \frac{{T_{{50}}}}{{(\ln 2)^{{1/SS}}}}
$$

M1/M2 estimate the scale $TD$ directly (code: `tvalpha`/`Alpha`); M3/M4 estimate $T_{{50}}$ (code: `tvt50`),
which is converted to $TD$. The shape $SS$ is `tvgamma`/`Gamma` in the code.

**Implementation note.** Long-acting injectable models commonly use Pumas `@delay`; the felodipine code writes the same input as a
depot with the Weibull hazard $h(\tau) = \frac{{SS}}{{TD}}(\tau/TD)^{{SS-1}}$:
$\frac{{dA_d}}{{dt}} = -h(\tau)A_d$, so $A_d(t) = \text{{Dose}}\cdot\exp[-(\tau/TD)^{{SS}}]$ and
$h(\tau)A_d(t) = \text{{Dose}}\cdot f_W(\tau) = I(t)$ exactly. In code, $\tau$ is floored at $10^{{-8}}$ h (M1, M2)
or $0.01$ h (M3, M4) because the density is singular at $\tau = 0$ when $SS < 1$ (clamping $SS \ge 1$ is not possible here: the felodipine arms span shapes 0.5–1.5).

**`@delay` form (M3d, M5, M6).** These models use Pumas `@delay` directly, as long-acting injectable models commonly do:
the dose record goes into the receiving compartment with bioavailability 0 and the input is
`@delay(Tlag + truncated(Weibull(SS, TD), τ0, Inf), cmt)` with $\tau_0 = 10^{{-4}}$ h, i.e.

$$
I(t) = \text{{Dose}}\cdot \frac{{f_W(\tau;\, SS, TD)}}{{1 - F_W(\tau_0;\, SS, TD)}}\quad (\tau \ge \tau_0),
\qquad 0 \text{{ otherwise}}.
$$

Three properties of `@delay` found while building it (prototypes 2026-10-24):

- `@dosecontrol` lags are ignored by `@delay`, so the lag has to sit inside the delay distribution.
- An untruncated `Tlag + Weibull`, or a dose record placed at $t = t_{{\text{{lag}}}}$, returns an Inf likelihood
  for $SS < 1$ (the density is singular at the dose time).
- The truncation at $\tau_0$ removes at most 0.4 % of the release (at $SS = 0.5$). It plays the part of the
  clock floor of the hazard form and is closer to the exact Weibull than the 0.01 h floor of M3.

With the lag fixed, M3d reproduces M3 on the 5 diagnostic arms (best LL within 3 units, identical Tmax and mean
curves). Estimating the lag in the `@delay` form fails (Inf at the start for $SS < 1$), so $t_{{\text{{lag}}}}$ stays
fixed in all models. Run times are similar to the hazard ODE (diagnostic fits: median 34–149 s per fit vs 39–127 s).

**Inter-individual variability.** $P_i = \theta_P\, e^{{\eta_{{P,i}}}}$, $\;\boldsymbol{{\eta}}_i \sim \mathcal{{N}}(\mathbf{{0}}, \Omega)$, $\Omega$ diagonal.

**Residual error.** $y_{{ij}} = C_{{ij}}\,(1 + \varepsilon_{{ij}})$, $\;\varepsilon_{{ij}} \sim \mathcal{{N}}(0, \sigma^2)$
(implemented as SD $= \sigma C + 10^{{-6}}$).

**Output.** $C = A_c / V_c$. All volumes and clearances are apparent (/F). FOCE, two starts per fit.

# M1 — 1-compartment + Weibull release (v1 baseline, 2026-10-16)

$$
\begin{{aligned}}
  \frac{{dA_c}}{{dt}} &= I(t) - \frac{{CL}}{{V_c}}\,A_c
\end{{aligned}}
$$

$\theta = (TD, SS, CL, V_c)$; IIV on $TD$, $CL$, $V_c$; 8 estimated parameters (4 $\theta$, 3 $\omega^2$, $\sigma$).

```julia
{M1}
```

# M2 — 2-compartment + Weibull release (2026-10-16, rejected)

$$
\begin{{aligned}}
  \frac{{dA_c}}{{dt}} &= I(t) - \frac{{CL + Q}}{{V_c}}\,A_c + \frac{{Q}}{{V_p}}\,A_p \\
  \frac{{dA_p}}{{dt}} &= \frac{{Q}}{{V_c}}\,A_c - \frac{{Q}}{{V_p}}\,A_p
\end{{aligned}}
$$

$\theta = (TD, SS, CL, V_c, Q, V_p)$; IIV on $TD$, $CL$, $V_c$; 10 estimated parameters.
On 2026-10-16, $Q$ or $V_p$ ran to their bounds from either start (peripheral compartment not identified).

```julia
{M2}
```

# M3 — Weibull release → first-order absorption → 2-compartment (v2, per arm)

$$
\begin{{aligned}}
  \frac{{dA_g}}{{dt}} &= I(t) - k_a\,A_g \\
  \frac{{dA_c}}{{dt}} &= k_a\,A_g - \frac{{CL + Q}}{{V_c}}\,A_c + \frac{{Q}}{{V_p}}\,A_p \\
  \frac{{dA_p}}{{dt}} &= \frac{{Q}}{{V_c}}\,A_c - \frac{{Q}}{{V_p}}\,A_p
\end{{aligned}}
$$

$$
T_{{50,i}} = \theta_{{T_{{50}}}}\, e^{{\eta_{{T_{{50}},i}}}}, \qquad
CL_i = \theta_{{CL}}\, e^{{\eta_{{CL,i}}}}, \qquad
V_{{c,i}} = \theta_{{V_c}}\, e^{{\eta_{{V_c,i}}}}
$$

$\theta = (T_{{50}}, SS, k_a, CL, V_c, Q, V_p)$; 11 estimated parameters (7 $\theta$, 3 $\omega^2$, $\sigma$).
Starts: $T_{{50}}$ and $SS$ at the arm's PK-Sim values, plus the two disposition starts below.
The gut compartment separates release from appearance in plasma, giving the PBPK profile's fast rise and plateau
without forcing release to be fast.

```julia
{M3}

{INITS}
```

# M4 — M3 paired: reference and test in the same subjects (prototype)

Reference dosed at $t = 0$, test at $t = 200$ h, each through its own Weibull input:
$I_R(t) = \text{{Dose}}\cdot f_W(t - t_{{\text{{lag}}}};\, SS_R, TD_R)$ and
$I_T(t) = \text{{Dose}}\cdot f_W(t - 200 - t_{{\text{{lag}}}};\, SS_T, TD_T)$ (code: `Depot1`, `Depot2`).

$$
\begin{{aligned}}
  \frac{{dA_g}}{{dt}} &= I_R(t) + I_T(t) - k_a\,A_g \\
  \frac{{dA_c}}{{dt}} &= k_a\,A_g - \frac{{CL + Q}}{{V_c}}\,A_c + \frac{{Q}}{{V_p}}\,A_p \\
  \frac{{dA_p}}{{dt}} &= \frac{{Q}}{{V_c}}\,A_c - \frac{{Q}}{{V_p}}\,A_p
\end{{aligned}}
$$

$TD_R = T_{{50,R}}/(\ln 2)^{{1/SS_R}}$, $TD_T = T_{{50,T}}/(\ln 2)^{{1/SS_T}}$. One subject-level $\eta_{{T_{{50}}}}$
scales both $T_{{50,R}}$ and $T_{{50,T}}$; $k_a$, $CL$, $V_c$, $Q$, $V_p$ and $\eta_{{CL}}$, $\eta_{{V_c}}$ are shared
across periods. 13 estimated parameters. Its AIC is not comparable with the per-arm models (two arms of data);
compare its LL with the sum of the M3 reference-arm and test-arm LLs.

```julia
{TP2}

{M4}
```

# M3d — M3 in the `@delay` form (2026-10-24)

The same equations as M3, with $I(t)$ written as the truncated `@delay` input above (dose record into `Gut`,
bioavailability 0). 11 estimated parameters, the M3 starts. Used to check that the delay-style input gives
the same fit as the hazard ODE before building M5 and M6 in that form.

```julia
{M3D}
```

# M5 — dual Weibull release → first-order absorption → 2-compartment (final model)

M3 peaks 0.4–0.7 h early on most arms and falls too fast at high shape: one release timescale and one $k_a$
cannot give both the fast initial phase and the prolonged late release of the PBPK profile. M5 replaces the
single Weibull with a two-Weibull mixture (`FF*dd1 + (1-FF)*dd2`):

$$
I(t) = \text{{Dose}}\left[FF\cdot \tilde f_W(\tau;\, SS_a, TD_a) + (1 - FF)\cdot \tilde f_W(\tau;\, SS_b, TD_b)\right],
\qquad TD_k = \frac{{T_{{50,k}}}}{{(\ln 2)^{{1/SS_k}}}}
$$

where $\tilde f_W$ is the Weibull density truncated at $\tau_0$ (see the `@delay` form above), and

$$
\begin{{aligned}}
  \frac{{dA_g}}{{dt}} &= I(t) - k_a\,A_g \\
  \frac{{dA_c}}{{dt}} &= k_a\,A_g - \frac{{CL + Q}}{{V_c}}\,A_c + \frac{{Q}}{{V_p}}\,A_p \\
  \frac{{dA_p}}{{dt}} &= \frac{{Q}}{{V_c}}\,A_c - \frac{{Q}}{{V_p}}\,A_p
\end{{aligned}}
$$

$$
T_{{50a,i}} = \theta_{{T_{{50a}}}}\, e^{{\eta_{{T_{{50}},i}}}}, \quad
T_{{50b,i}} = \theta_{{T_{{50b}}}}\, e^{{\eta_{{T_{{50}},i}}}}, \quad
CL_i = \theta_{{CL}}\, e^{{\eta_{{CL,i}}}}, \quad
V_{{c,i}} = \theta_{{V_c}}\, e^{{\eta_{{V_c,i}}}}
$$

$\theta = (T_{{50a}}, SS_a, T_{{50b}}, SS_b, FF, k_a, CL, V_c, Q, V_p)$, $0 \le FF \le 1$; one release $\eta$ scales
both components; $t_{{\text{{lag}}}} = 0.3113$ h fixed. 14 estimated parameters (10 $\theta$, 3 $\omega^2$, $\sigma$).

**Starts.** Two starts bracket the arm's PK-Sim $T_{{50}}$ ($T$) and shape ($s$). The 27-arm production fits
(`02c_fit_poppk_m5.jl`) use the absorption-limited start only: on the 7 arms fitted from both starts it was never
worse than the release-limited one and had the higher LL on 3. Fits take about 30 min per arm single-threaded.

```julia
{M5STARTS}
```

**Interpretation.** The two components are a flexible description of the release curve, not two physical
populations: $FF$ ranges 0.02–0.97 across arms, and on some arms component $a$ is a small fast burst
($T_{{50a}}$ 0.2–0.4 h). For comparison with the PK-Sim input, the per-arm parameter table reports summaries of
the mixture: `tvt50` = time of 50 % released (lag excluded) and `tvgamma` = shape of the single Weibull through the
mixture's 25 % and 75 % release times. Like M3's $T_{{50}}$, the summary $T_{{50}}$ sits well above the PK-Sim
dissolution $T_{{50}}$ (≈ 8–9 h at PK-Sim 4–5 h) because the PopPK input lumps release, gut transit and
absorption. At the reference arm $FF$ (RSE 45 %) and $Q$, $V_p$ (RSE 55–70 %) are poorly determined; $V_c$
reaches ≈ 2300 L at shape 1.5, as in M3.

```julia
{M5}
```

# M6 — Weibull release → 3 transit compartments → 2-compartment (tested, not selected)

$$
\begin{{aligned}}
  \frac{{dA_{{1}}}}{{dt}} &= I(t) - k_{{tr}}\,A_{{1}}, \qquad
  \frac{{dA_{{k}}}}{{dt}} = k_{{tr}}\,(A_{{k-1}} - A_{{k}}), \; k = 2, 3 \\
  \frac{{dA_c}}{{dt}} &= k_{{tr}}\,A_{{3}} - \frac{{CL + Q}}{{V_c}}\,A_c + \frac{{Q}}{{V_p}}\,A_p, \qquad
  \frac{{dA_p}}{{dt}} = \frac{{Q}}{{V_c}}\,A_c - \frac{{Q}}{{V_p}}\,A_p
\end{{aligned}}
$$

$I(t)$ as M3d; the transit chain replaces $k_a$. 11 estimated parameters. It fixes the shape-1.5 curve
(LL 19 vs M3 −59) but is no better than M3 on the other diagnostic arms, and 2 of its 10 fits failed.

```julia
{M6}
```

# Model selection

**Diagnostic arms** (4.78/0.764 reference, 4.78/0.5, 4.78/1.5, 2.39/0.764, 9.56/0.764; best start per arm,
mean curves and per-subject Tmax against the noise-free PK-Sim profiles):

| Model | Best LL by arm (ref, 0.5, 1.5, 2.39, 9.56) | Tmax bias (median per subject) | Mean-curve error 0.5–12 h |
|---|---|---|---|
| M3 | −354, −487, −59, −481, −108 | −0.7 to +0.1 h | 2.0–6.6 % |
| M3d | −354, −488, −56, −481, −107 | as M3 | as M3 |
| M3, no lag (`@delay`) | −755, −754, −48, −866, −500 | +60 % overshoot at 0.5 h | 2.9–8.4 % |
| M3, lag estimated (ODE) | −526, −762, −329, −808, −106 | lag ≈ 0.27–0.33 h, LL worse | — |
| **M5** | **−348, −439, +38, −433, −85** | **−0.3 to +0.15 h** | **1.2–4.6 %** |
| M6 | −365, −479, +19, −484, −109 | −0.45 to 0 h | 2.1–4.3 % |

**All 27 arms, M5 vs M3** (`felo_09_m5_vs_m3_all27.jl`; Cmax/AUC against the PBPK values in the same subjects,
T/R GMR against the 4.78/0.764 reference arm):

| | M3 (fallback) | M5 (final) |
|---|---|---|
| AIC better | 0 / 27 arms | 27 / 27 (median ΔAIC −48, range −289 to −6.5) |
| Cmax ratio PopPK/PBPK | 0.976–1.043 | 0.976–1.019 |
| AUC₀₋ₜ ratio PopPK/PBPK | 0.985–1.012 | 0.989–1.004 |
| Cmax T/R GMR error vs PBPK | max 3.5 %, mean 1.6 % | max 3.3 %, mean 1.3 % |
| Tmax bias vs PK-Sim (median of arms) | 0.45 h | 0.20 h |
| Subjects with Tmax within 0.5 h | 64 % | 76 % |
| 2-D safe-space T50 limits vs PBPK | mean \|Δ\| 0.08 h, max 0.20 h, bias +0.03 h | mean \|Δ\| 0.09 h, max 0.25 h, bias −0.09 h |

M5 describes each arm better, but its safe-space limits are not closer to PBPK than M3's. M5 overpredicts the
reference-arm Cmax by 0.9 % while its test arms are nearly unbiased, so its T/R GMRs run about 1.5 % low at shape
≥ 0.9 and the limits shift left. M3 overpredicts Cmax on the reference and test arms alike (up to 3.5 % at shape 1.2)
and the two largely cancel. Either model leaves the BE conclusions unchanged (all limits within 0.25 h of PBPK).
One M5 arm (6.31/1.2) still peaks early (median −0.85 h).

**BE decision vs PBPK** (`felo_04_pbpk_poppk_compare.jl`; 1-D axes through the reference, WSV fasted,
Type I from 20 000 simulated trials, seed 1234). Cmax sets every limit; AUC₀₋ₜ never reaches the BE bounds.

| | T50 axis lower limit (h) | T50 axis upper limit (h) | Type I at declared limit, Cmax, n = 24 / 54 (lower; upper) |
|---|---|---|---|
| PBPK | 3.25 | 6.96 | 4.9 / 5.1 %; 4.7 / 5.0 % (nominal) |
| M3 | 3.40 | 6.83 | 7.6 / 10.2 %; 5.8 / 6.8 % |
| M5 | 3.27 | 6.86 | 5.3 / 5.8 %; 5.6 / 6.3 % |

On the shape axis (T50 4.78 h) PBPK loses BE at shape 1.465. Both PopPK models underpredict the shape-1.5 Cmax
GMR (M3 1.235, M5 1.227, PBPK 1.269), so neither reaches 1.25 within the tested range.
"""

(HERE / "felodipine_poppk_models.qmd").write_text(doc)
print("wrote", HERE / "felodipine_poppk_models.qmd")
