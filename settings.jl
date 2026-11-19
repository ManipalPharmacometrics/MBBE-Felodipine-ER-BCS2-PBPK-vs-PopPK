# =============================================================================
# settings.jl — every number the analysis uses, in one place.
# Each step script starts with include("settings.jl").
# =============================================================================

const HERE    = @__DIR__
const DATA    = joinpath(HERE, "data", "pksim_exports")
const RESULTS = joinpath(HERE, "results")

# ---- the drug and the study -------------------------------------------------
const MW        = 384.26                         # felodipine, g/mol
const DOSE_NMOL = 10.0 / MW * 1e6                # 10 mg ER tablet = 26,024 nmol

# Rich sampling schedule (h), 18 samples to 96 h. The dose is at 0 h.
const SCHEDULE = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 24, 36, 48, 72, 96]
const PROP_CV  = 0.15                            # 15% proportional assay error
const SEED_DATASET = 20261016                    # random seed for the assay error

# ---- formulations -----------------------------------------------------------
# PK-Sim simulated 27 ER formulations that differ in the Weibull release curve:
# T50 (h, time to 50% released) and shape. Reference = the marketed tablet.
const REF = (T50 = 4.78, shape = 0.764)

# ---- bioequivalence ---------------------------------------------------------
const BE_LOWER = 0.80
const BE_UPPER = 1.25
const ENDPOINTS = ["Cmax", "AUC_tEnd"]            # AUC_tEnd = AUC 0–96 h

# Within-subject CV, fasted 10 mg 2x2 crossover (UK MHRA public assessment
# report, Felodipine/Delofine XL, PL 20117/0397-0399): Cmax 35.7%, AUC0-t 25.3%.
const WSV_CV = Dict("Cmax" => 0.3571, "AUC_tEnd" => 0.2525)
# Reference within-subject CV for the scaled methods (FDA RSABE, EMA ABEL):
# ~34%, back-calculated from the ABEL limits of the replicate study (approximate).
const CV_WR = 0.34

const N_ABE   = 54          # 2x2 study powered at 80% for CV 35.7%, GMR 0.95
const N_SUBJ  = [24, 54]    # study sizes for the axis analyses (steps 5–6)
const N_SCALED = 24         # replicate study size for RSABE / ABEL
const SEED_TRIALS = 1234    # random seed for the virtual BE trials

const N_TRIAL_POWER = 1_000                      # trials per formulation / grid point
const N_TRIAL_TYPE1 = 20_000                     # SE of a 5% rate ~0.15 percentage points
const POWER_TARGET  = 0.80
const N_TRIAL_CHECK = 2_000                      # step 6b: trials per point run through both tests

# The two 1-D slices of the design used in steps 5–6
#   T50 axis  : shape fixed at 0.764, 9 formulations, T50 2.39–9.56 h
#   shape axis: T50 fixed at 4.78 h, 7 formulations, shape 0.5–1.5
# Between formulations, values are interpolated in log(value) against log(T50)
# on the T50 axis and against shape on the shape axis. 161-point grids.
const AXES = [(name = "T50",   col = :T50,   fixed = :shape, fixed_value = REF.shape, ref = REF.T50,   log_x = true),
              (name = "shape", col = :shape, fixed = :T50,   fixed_value = REF.T50,   ref = REF.shape, log_x = false)]
const GRID_POINTS = 161
