# Model-Based Bioequivalence: Felodipine ER (BCS II), PBPK vs PopPK

Do PBPK and PopPK reach the same bioequivalence (BE) decisions for a dissolution-limited (BCS II)
extended-release product? This repository extends the metformin XR comparison to felodipine ER 10 mg,
sweeping release over two quality attributes at once: the time to 50% release (T50) and the shape of
the Weibull release curve.

- **PBPK (the reference).** The Open Systems Pharmacology felodipine model, run in PK-Sim 12.2, with
  27 ER formulations on a T50 × shape grid, each simulated in a virtual population.
- **PopPK (the method under test).** Pumas models fitted per formulation to rich simulated studies
  (18 samples per subject, 15% proportional error). Six structures were tested (M1–M6). **M5 is the
  final model** (two Weibull release components, first-order absorption, two-compartment disposition);
  M3 is the fallback.
- **Output.** The BE safe space on the T50 × shape plane from each method, type I error at each
  method's own boundary, power, and RSABE/ABEL sensitivity.

## Summary of results

Both methods are calibrated at their own boundary (type I error 4.5–5.1% for Cmax and the joint
endpoint at n = 54). PopPK, however, moves the GMR boundaries: on the T50 axis the fast-side boundary
shifts by −15% and the slow-side boundary by −7%, against −1.1% and −6.8% for metformin. Agreement
between PopPK and PBPK BE decisions therefore needs an adequate structural model, and for felodipine
it fails away from the reference formulation even with rich sampling.

## Layout

```
docs/                       start here: full work summary (pptx/html), PK-Sim work summary (docx),
                            discussion and explainer decks (offline html), model testing log,
                            PopPK model equations (pdf)
Felodipine-BCSII/           main project folder
  reference/                OSP model snapshot, simulation overrides, digitised literature data
  pksim-particle-size-pk/   particle-size scan of the IR tablet (attribute later dropped)
  pksim-releaserate-cqa-er/ PK-Sim exports: model checks, colonic CYP3A4 scan, the 27-arm sweep
                            (matrix_Sweep_td50_shape/), population runs, model snapshots
  scripts/                  felo_02 … felo_09: BE sweep, safe space, PBPK vs PopPK, diagnostics,
                            model comparison (Julia); deck builder and spreadsheet readers
  output_sweep_*/           PBPK BE results (coarse, 2D, fed-CV sensitivity, 1994 tablets)
  poppk/                    PopPK results collected from the fits: output_rich (M1),
                            output_rich_v2 (M3), output_rich_m5 (M5), sweeps through the same
                            BE analysis, model_comparison/ (tables, equations, testing log)
  slides/                   deck sources
poppk-fits/                 Pumas models (src/), fit scripts (scripts/) and the simulated dataset
julia-env/                  Julia environment the PopPK fits ran in
programs/                   Julia environment for the BE analysis, plus julia/be_utils.jl
```

## PK-Sim model

The PBPK model is the OSP library felodipine model
([Open-Systems-Pharmacology/Felodipine-Model](https://github.com/Open-Systems-Pharmacology/Felodipine-Model),
master snapshot of 14 Oct 2026; Bergmann et al., *Pharmaceutics* 2022;14:1474). The ER project used
here is stored as a snapshot:
`Felodipine-BCSII/pksim-releaserate-cqa-er/Felodipine-ER-colon50-15102026-snapshot.json`
(in PK-Sim, File → Load from snapshot).

When rebuilding, use **Clone**, not **Configure**: Configure drops the OSP simulation-level overrides
(22 mucosal permeability values of 0.0891674187 cm/min and the tablet time-delay factor 0.0695415108).
The full list is in `Felodipine-BCSII/reference/ER_tablet_10mg_sim_overrides.csv`.

PK-Sim is not needed to rerun the analysis: all 27 formulations are already exported.

## Running it

Julia 1.11 was used. Scripts build their paths from their own folder, so the repository can be
cloned anywhere.

**BE analysis (open-source Julia packages only)**, from `programs/`:

```bash
cd programs
julia --project=. -e 'using Pkg; Pkg.instantiate()'
F=../Felodipine-BCSII
julia --project=. $F/scripts/felo_02_be_sweep.jl $F/pksim-releaserate-cqa-er/matrix_Sweep_td50_shape $F/output_sweep_2d
julia --project=. $F/scripts/felo_03_safespace_2d.jl $F/output_sweep_2d
julia --project=. $F/scripts/felo_02_be_sweep.jl $F/poppk/output_rich_m5/pk_analyses_format $F/poppk/output_sweep_2d_poppk_m5
julia --project=. $F/scripts/felo_03_safespace_2d.jl $F/output_sweep_2d $F/poppk/output_sweep_2d_poppk_m5 PopPK-M5
julia --project=. $F/scripts/felo_04_pbpk_poppk_compare.jl $F/poppk/output_rich_m5
julia --project=. $F/scripts/felo_05_poppk_diagnostics.jl $F/poppk/output_rich_m5
```

These commands overwrite the stored outputs; give `felo_02`/`felo_03` a new output folder to compare.
`felo_04` runs 20,000 type I trials per point; set `N_T1=500 N_PW=100` for a quick check. Seeds are
fixed, so a full run reproduces the stored numbers.

**PopPK fits (needs a licensed Pumas installation)**, from `poppk-fits/`:

```bash
cd poppk-fits
J="julia -J$PUMAS_SO --project=../julia-env"
$J scripts/01_build_dataset.jl
POPPK_OUT=$PWD/output_m5 $J scripts/02c_fit_poppk_m5.jl
POPPK_OUT=$PWD/output_m5 LEARN_DST=$PWD/../Felodipine-BCSII/poppk/output_rich_m5_new $J scripts/03_collect.jl
```

`julia-env/Project.toml` lists more packages than the felodipine scripts need. If they do not resolve,
Pumas plus CSV, DataFrames, Distributions and CairoMakie is enough.

## Not included

- The raw per-arm Pumas fit outputs (about 180 MB, including serialized fits). The collected results
  are in `Felodipine-BCSII/poppk/`; the raw outputs are regenerated by the fit scripts above.
  `felo_07` and `felo_08` read the raw outputs, so run the fits before those two.
- Individual predicted profiles on a fine time grid (`poppk_ipred_thin.csv`, 18 MB each), written by
  the fit scripts.
- Published papers and product monographs used as references; they are cited in
  `Felodipine-BCSII/reference/*.md`.

All data in this repository are simulated, or digitised from published figures. No clinical or
patient-level data are included.

## Related work

- [MBBE-Metformin-BCS3-PBPK-vs-PopPK](https://github.com/ManipalPharmacometrics/MBBE-Metformin-BCS3-PBPK-vs-PopPK):
  the first product in this comparison.
- [Uncertainty Propagation in MBBE](https://github.com/ManipalPharmacometrics/Uncertainity-Propogation-in-MBBE).

## Author

Nikita Ramwani
