# =============================================================================
# be_utils.jl
# Shared helpers for the Metformin XR PBPK virtual-bioequivalence (VBE) and
# safe-space analysis. PRIMARY implementation (Julia / Pumas ecosystem).
#
# Project: PBPK vs PopPK for BE safe space + Type I error. Reference formulation: T50 = 3.0 h (as specified).
#
# Design note -- PAIRED crossover:
#   PK-Sim simulates the SAME virtual individuals (IndividualId 0..N-1) under
#   every formulation, so test-vs-reference is a fully paired 2x2 crossover with
#   no period/sequence effect. We analyse within-subject log-differences, the
#   textbook average-BE (TOST) estimator.
#
# BE package note:
#   Pumas `Bioequivalence`/`NCA` (licensed, in PumasRegistry) can be dropped in
#   once the env is instantiated with a JuliaHub token. The open-registry
#   `Bioequivalence` UUID is broken (does not define its module), so the TOST is
#   implemented here directly. The R scripts under programs/R provide an
#   independent cross-check of these numbers.
#
# Multi-drug note (2026-10-15, felodipine ER = MBBE drug two):
#   The metformin API (REF_T50, WSV_CV, parse_t50, load_all_formulations(dir))
#   is unchanged. Other drugs pass a `drug` key: DRUG_CONFIG / ref_formulation,
#   wsv_cv(drug, arm, endpoint), load_all_formulations(dir; drug) with (T50,
#   shape) parsed from filenames, sample_size_abe_2x2, and full-replicate
#   reference scaling (FDA RSABE / EMA ABEL) via vbe_replicate_scaled.
# =============================================================================

module BEUtils

using CSV, DataFrames, Distributions, Random, Statistics

export REF_T50, BE_LOWER, BE_UPPER, BE_CI_LEVEL, BE_METRICS,
       WSV_CV, WSV_CV_IR, N_TRIAL_DEFAULT, N_SUBJ_LIST,
       load_all_formulations, be_paired, be_all_metrics, vbe_replicate,
       vbe_replicate_joint,
       # multi-drug (felodipine) extensions -- see "Drug configuration" below
       DRUG_CONFIG, WSV_TABLE, wsv_cv, ref_formulation, parse_formulation,
       reference_rows, power_abe_2x2, sample_size_abe_2x2,
       scaled_be_replicate, vbe_replicate_scaled

# ---- Configuration ----------------------------------------------------------
const REF_T50      = 3.0          # reference formulation T50 (h)
const BE_LOWER     = 0.80         # FDA/EMA average-BE limits (log scale)
const BE_UPPER     = 1.25
const BE_CI_LEVEL  = 0.90         # 90% CI <=> two one-sided tests at alpha=0.05

# PK metric label => PK-Sim `Parameter` name.
#   AUC_tEnd (AUC_0-t) is preferred for XR / high-extrapolation products;
#   AUC_inf reported alongside.
const BE_METRICS = ["Cmax" => "C_max", "AUC_inf" => "AUC_inf",
                    "AUC_tEnd" => "AUC_tEnd"]

# Endpoint-specific within-subject (intra-subject) CVs for the VBE simulation.
# Metformin extended-release, sourced:
#   FDA NDA 212595 ClinPharm review: Cmax 12.2%, AUC_0-t 11.2%, AUC_inf 11.4%
#   (corroborated Cmax 14.3% / AUC 10.3%, Chung et al. 2018 TCP 26(1):6-9).
const WSV_CV = Dict("Cmax" => 0.122, "AUC_inf" => 0.114, "AUC_tEnd" => 0.112)

# High-variability sensitivity scenario — metformin IMMEDIATE-release / upper
# bound: Cmax 23.7% (IR tablets, lit. range 22.7-25.7%; observed up to 24.8%),
# AUC 14.7% (pooled 142-study mean; metformin IR study reported 13.8%).
const WSV_CV_IR = Dict("Cmax" => 0.237, "AUC_inf" => 0.147, "AUC_tEnd" => 0.147)

# VBE simulation defaults (as specified): 500 replicate trials, study sizes
# of 48 and 86 subjects. NB the PBPK population has 24 unique virtual subjects,
# so n>24 trials bootstrap-with-replacement from those 24 physiologies.
const N_TRIAL_DEFAULT = 500
const N_SUBJ_LIST     = [48, 86]

# ---- Drug configuration (multi-drug MBBE) -------------------------------------
# The metformin constants above are kept verbatim so the metformin scripts and
# their outputs are unchanged. New analyses pass a `drug` key instead.
#
# ref_T50 / ref_shape : reference Weibull release (h, unitless). For metformin
#   the shape is fixed and absent from filenames, so it stays `missing`.
# Felodipine reference (Felodipine-BCSII project notes §9.6, §9.14):
#   shipped OSP ER Weibull T50 286.95 min = 4.78 h, shape 0.764, lag 18.68 min,
#   with colonic CYP3A4 at 50% of the small-intestine level.
const DRUG_CONFIG = Dict(
    :metformin  => (ref_T50 = REF_T50, ref_shape = missing, default_arm = :fasted),
    :felodipine => (ref_T50 = 4.78,    ref_shape = 0.764,   default_arm = :fasted),
)

# Within-subject CV by drug x arm x endpoint.
# Felodipine (scope planning notes §4.3; UK MHRA PAR
# Felodipine/Delofine XL PL 20117/0397-0399, 10 mg 2x2 crossovers):
#   fasted Cmax 35.71%, AUC0-t 25.25%, AUC0-inf 24.34%   (Study 1)
#   fed    Cmax 36.51%, AUC0-t 21.92%, AUC0-inf 22.16%   (Study 2)
#   CVwR for the Cmax RSABE sensitivity ~ 0.34 (Study 4, back-calculated from
#   the EMA ABEL limits 77.95-128.28%; approximate).
const WSV_TABLE = Dict(
    (:metformin,  :fasted)          => WSV_CV,
    (:metformin,  :ir_sensitivity)  => WSV_CV_IR,
    (:felodipine, :fasted) => Dict("Cmax" => 0.3571, "AUC_inf" => 0.2434, "AUC_tEnd" => 0.2525),
    (:felodipine, :fed)    => Dict("Cmax" => 0.3651, "AUC_inf" => 0.2216, "AUC_tEnd" => 0.2192),
    (:felodipine, :cvwr_rsabe)      => Dict("Cmax" => 0.34),
)

"Within-subject CV map (endpoint => CV) for `drug` and `arm`."
function wsv_cv(drug::Symbol, arm::Symbol = DRUG_CONFIG[drug].default_arm)
    haskey(WSV_TABLE, (drug, arm)) || error("No WSV_CV for ($drug, $arm)")
    WSV_TABLE[(drug, arm)]
end
wsv_cv(drug::Symbol, arm::Symbol, endpoint::AbstractString) = wsv_cv(drug, arm)[endpoint]

"Reference formulation `(T50 = h, shape = β)` for `drug`."
ref_formulation(drug::Symbol) =
    (T50 = DRUG_CONFIG[drug].ref_T50, shape = DRUG_CONFIG[drug].ref_shape)

# ---- Data loading ------------------------------------------------------------

"Extract T50 (h) from a PK-Sim filename like `..._T50_3.5h_Simpop-...csv`."
function parse_t50(path::AbstractString)
    m = match(r"T50_([0-9.]+)h", basename(path))
    m === nothing && error("No T50 token in $(basename(path))")
    parse(Float64, m.captures[1])
end

"""
    read_pk_analyses(path) -> DataFrame

Read one PK-Sim `*-PK-Analyses.csv` (long: IndividualId, QuantityPath,
Parameter, Value, Unit) and return a wide per-individual table with one column
per metric in `BE_METRICS`, tagged with the formulation `T50`.
"""
function read_pk_analyses(path::AbstractString)
    raw = CSV.read(path, DataFrame; normalizenames = true)
    # normalizenames turns "IndividualId" etc. into symbols; Parameter has no spaces.
    long = filter(r -> r.Parameter in last.(BE_METRICS), raw)
    ids  = sort(unique(long.IndividualId))
    out  = DataFrame(T50 = fill(parse_t50(path), length(ids)), IndividualId = ids)
    for (label, pcol) in BE_METRICS
        sub = filter(r -> r.Parameter == pcol, long)
        d   = Dict(sub.IndividualId .=> Float64.(sub.Value))
        out[!, label] = [get(d, id, missing) for id in ids]
    end
    out
end

"""
    parse_formulation(path, drug) -> (T50 = h, shape = β)

Read the Weibull release parameters from a PK-Sim export name. Either token may
be absent and then defaults to the drug's reference value, so all of these work:
  `Metformin_T50_3.5h_Simpop-PK-Analyses.csv`  -> (3.5,  missing)
  `Pop_T50_2.39h-PK-Analyses.csv` (felodipine) -> (2.39, 0.764)
  `Pop_b0.5-PK-Analyses.csv`      (felodipine) -> (4.78, 0.5)
  `Pop_T50_3.35h_b0.9-PK-Analyses.csv`         -> (3.35, 0.9)
The shape token must carry a decimal point (`b1.0`, not `b10`) so that
`b0.1` and `b1.0` cannot be confused.
"""
function parse_formulation(path::AbstractString, drug::Symbol)
    name = basename(path)
    ref  = ref_formulation(drug)
    mt   = match(r"T50_([0-9]+(?:\.[0-9]+)?)h", name)
    mb   = match(r"_b([0-9]*\.[0-9]+)(?=[-_]|$)", name)
    occursin(r"_b[0-9]+(?=[-_]|$)", name) && mb === nothing &&
        error("Shape token without a decimal point in $name (write e.g. b0.764)")
    t50  = mt === nothing ? ref.T50   : parse(Float64, mt.captures[1])
    sh   = mb === nothing ? ref.shape : parse(Float64, mb.captures[1])
    (T50 = t50, shape = sh)
end

"""
    read_pk_analyses(path, drug; quantity = "PeripheralVenousBlood") -> DataFrame

Drug-aware reader. Felodipine exports carry several outputs (peripheral venous
plasma, venous blood), so rows are restricted to the `QuantityPath` containing
`quantity`. Adds `drug` and `shape` columns next to `T50`.
"""
function read_pk_analyses(path::AbstractString, drug::Symbol;
                          quantity = "PeripheralVenousBlood")
    raw  = CSV.read(path, DataFrame; normalizenames = true)
    long = filter(r -> r.Parameter in last.(BE_METRICS) &&
                       occursin(quantity, r.QuantityPath), raw)
    isempty(long) && error("No $quantity PK rows in $(basename(path))")
    f    = parse_formulation(path, drug)
    ids  = sort(unique(long.IndividualId))
    out  = DataFrame(drug = fill(drug, length(ids)), T50 = fill(f.T50, length(ids)),
                     shape = fill(f.shape, length(ids)), IndividualId = ids)
    for (label, pcol) in BE_METRICS
        sub = filter(r -> r.Parameter == pcol, long)
        nrow(sub) == length(ids) || error("Duplicate $pcol rows in $(basename(path))")
        d   = Dict(sub.IndividualId .=> Float64.(sub.Value))
        out[!, label] = [get(d, id, missing) for id in ids]
    end
    out
end

"""
    load_all_formulations(dir; drug = nothing)

Read & stack every `*-PK-Analyses.csv` under `dir`. Without `drug` this is the
original metformin loader (T50 only). With `drug`, rows also carry `drug` and
`shape`, and two files mapping to the same (T50, shape) raise an error.
"""
function load_all_formulations(dir::AbstractString; drug = nothing)
    files = filter(f -> endswith(f, "-PK-Analyses.csv"), readdir(dir; join = true))
    isempty(files) && error("No *-PK-Analyses.csv under $dir")
    drug === nothing && return reduce(vcat, read_pk_analyses.(files))
    keys = [parse_formulation(f, drug) for f in files]
    for k in unique(keys)
        dup = files[keys .== Ref(k)]
        length(dup) > 1 && error("Same formulation $k in: " * join(basename.(dup), ", "))
    end
    reduce(vcat, [read_pk_analyses(f, drug) for f in files])
end

"Rows of `pk` belonging to `drug`'s reference formulation (errors if absent)."
function reference_rows(pk::DataFrame, drug::Symbol)
    ref = ref_formulation(drug)
    sel = ismissing(ref.shape) ? (pk.T50 .== ref.T50) :
                                 (pk.T50 .== ref.T50) .& isequal.(pk.shape, ref.shape)
    any(sel) || error("Reference $ref not found for $drug")
    pk[sel, :]
end

# ---- Bioequivalence statistics ----------------------------------------------

"""
    be_paired(test, ref; ...) -> NamedTuple

Paired (crossover) average-BE for ONE metric: within-subject log-differences,
GMR, 90% CI (t-based, TOST), and BE pass/fail. Subjects with non-finite or
non-positive values (e.g. AUC_inf extrapolation failures) are dropped.
"""
function be_paired(test, ref; ci_level = BE_CI_LEVEL,
                   lower = BE_LOWER, upper = BE_UPPER)
    ok = @. !ismissing(test) & !ismissing(ref) &
            isfinite(test) & isfinite(ref) & (test > 0) & (ref > 0)
    d  = log.(test[ok]) .- log.(ref[ok])
    n  = length(d)
    md = mean(d)
    se = std(d) / sqrt(n)
    tq = quantile(TDist(n - 1), 1 - (1 - ci_level) / 2)
    lo = exp(md - tq * se)
    hi = exp(md + tq * se)
    (n = n, GMR = exp(md), CI_lower = lo, CI_upper = hi,
     BE = (lo >= lower) && (hi <= upper))
end

"Run `be_paired` for every metric of one test formulation vs the reference."
function be_all_metrics(test_df::DataFrame, ref_df::DataFrame)
    t = sort(test_df, :IndividualId)
    r = sort(ref_df,  :IndividualId)
    @assert t.IndividualId == r.IndividualId "subject IDs must align (paired design)"
    rows = map(first.(BE_METRICS)) do m
        res = be_paired(t[!, m], r[!, m])
        merge((metric = m,), res)
    end
    DataFrame(rows)
end

# ---- Virtual BE replicate trials (power / Type I error) ----------------------

"""
    vbe_replicate(test, ref; n_subj, n_trial, wsv_cv, seed) -> NamedTuple

Resample `n_subj`-subject crossover trials (bootstrap, with optional
within-subject residual CV `wsv_cv`) and return the BE pass fraction:
  * test != reference -> POWER (probability of correctly concluding BE)
  * test == reference -> TYPE I ERROR, but only meaningful with wsv_cv > 0
    (identical paired PBPK values give a degenerate self-comparison).
`wsv_cv` is an ASSUMPTION; set it from metformin within-subject CV literature.
"""
function vbe_replicate(test, ref; n_subj = nothing, n_trial = 2000,
                       wsv_cv = 0.0, seed = 1234,
                       ci_level = BE_CI_LEVEL, lower = BE_LOWER, upper = BE_UPPER)
    rng  = MersenneTwister(seed)
    keep = @. !ismissing(test) & !ismissing(ref) &
              isfinite(test) & isfinite(ref) & (test > 0) & (ref > 0)
    tv0  = Float64.(test[keep]); rv0 = Float64.(ref[keep])
    N    = length(rv0)
    nss  = n_subj === nothing ? N : n_subj
    sw   = sqrt(log(wsv_cv^2 + 1))          # log-scale SD from CV (0 -> 0)
    npass = 0; gmrs = Vector{Float64}(undef, n_trial)
    for t in 1:n_trial
        idx = rand(rng, 1:N, nss)
        tv  = tv0[idx]; rv = rv0[idx]
        if sw > 0
            tv = tv .* exp.(randn(rng, nss) .* sw)
            rv = rv .* exp.(randn(rng, nss) .* sw)
        end
        res = be_paired(tv, rv; ci_level, lower, upper)
        npass += res.BE
        gmrs[t] = res.GMR
    end
    (n_trial = n_trial, n_subj = nss, wsv_cv = wsv_cv,
     pass_rate = npass / n_trial, GMR_med = median(gmrs))
end

"""
    vbe_replicate_joint(testC, refC, testA, refA; cv_C, cv_A, ...) -> NamedTuple

Joint VBE: a trial concludes BE only if BOTH Cmax and AUC pass. Returns the
joint pass-rate. `cv_C`, `cv_A` are the metric-specific within-subject CVs.

Within-subject CV defaults (metformin ER, intra-subject):
  Cmax = 0.122, AUC0-t = 0.112  — FDA NDA 212595 ClinPharm review
  (corroborated: Cmax 0.143 / AUC 0.103, Chung et al. 2018 TCP 26(1):6-9).
"""
function vbe_replicate_joint(testC, refC, testA, refA;
                             n_subj = nothing, n_trial = 2000,
                             cv_C = 0.122, cv_A = 0.112, seed = 1234,
                             ci_level = BE_CI_LEVEL, lower = BE_LOWER, upper = BE_UPPER)
    rng  = MersenneTwister(seed)
    keep = @. !ismissing(testC) & !ismissing(refC) & !ismissing(testA) & !ismissing(refA) &
              isfinite(testC) & isfinite(refC) & isfinite(testA) & isfinite(refA) &
              (testC > 0) & (refC > 0) & (testA > 0) & (refA > 0)
    tc = Float64.(testC[keep]); rc = Float64.(refC[keep])
    ta = Float64.(testA[keep]); ra = Float64.(refA[keep])
    N   = length(tc)
    nss = n_subj === nothing ? N : n_subj
    swC = sqrt(log(cv_C^2 + 1)); swA = sqrt(log(cv_A^2 + 1))
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, nss)
        tcv = tc[idx]; rcv = rc[idx]; tav = ta[idx]; rav = ra[idx]
        if swC > 0
            tcv = tcv .* exp.(randn(rng, nss) .* swC); rcv = rcv .* exp.(randn(rng, nss) .* swC)
        end
        if swA > 0
            tav = tav .* exp.(randn(rng, nss) .* swA); rav = rav .* exp.(randn(rng, nss) .* swA)
        end
        rC = be_paired(tcv, rcv; ci_level, lower, upper)
        rA = be_paired(tav, rav; ci_level, lower, upper)
        npass += (rC.BE && rA.BE)
    end
    (n_trial = n_trial, n_subj = nss, cv_C = cv_C, cv_A = cv_A,
     pass_rate = npass / n_trial)
end

# ---- Sample size (2x2 crossover, average BE) ---------------------------------

"""
    power_abe_2x2(cv, gmr, n; lower, upper, alpha = 0.05) -> Float64

TOST power for a balanced 2x2 crossover with `n` subjects in total, using the
noncentral-t approximation (PowerTOST `method = "nct"`): SE = s_w * sqrt(2/n),
df = n - 2, s_w = sqrt(log(cv^2 + 1)).
"""
function power_abe_2x2(cv, gmr, n; lower = BE_LOWER, upper = BE_UPPER, alpha = 0.05)
    df = n - 2
    se = sqrt(log(cv^2 + 1)) * sqrt(2 / n)
    tq = quantile(TDist(df), 1 - alpha)
    p  = cdf(NoncentralT(df, (log(gmr) - log(upper)) / se), -tq) -
         cdf(NoncentralT(df, (log(gmr) - log(lower)) / se),  tq)
    max(p, 0.0)
end

"Smallest even total `n` (>= 12) reaching `power` for a 2x2 ABE study."
function sample_size_abe_2x2(cv; gmr = 0.95, power = 0.80, kwargs...)
    n = 12
    while power_abe_2x2(cv, gmr, n; kwargs...) < power
        n += 2
        n > 10_000 && error("No n <= 10000 reaches power $power")
    end
    n
end

# ---- Reference-scaled BE (Cmax sensitivity for highly variable drugs) --------

const RSABE_SIGMA_W0 = 0.25                        # FDA regulatory constant
const RSABE_THETA    = (log(1.25) / RSABE_SIGMA_W0)^2
const RSABE_SWR_MIN  = 0.294                       # FDA: scale only if sWR >= 0.294
const ABEL_K         = 0.760                       # EMA widening constant
const ABEL_SWR_CAP   = sqrt(log(0.50^2 + 1))       # EMA: widening capped at CVwR 50%

"""
    scaled_be_replicate(t1, t2, r1, r2; method = :fda_rsabe, alpha = 0.05) -> NamedTuple

One full-replicate (TRTR/RTRT) trial on log-scale per-subject values. There are
no period or sequence effects in the virtual data, so the subject-level
estimators are used: I = mean(T) - mean(R), D = R1 - R2, sWR^2 = var(D)/2,
df = n - 2 (two sequences, as in the FDA guidance).

`:fda_rsabe`: if sWR >= 0.294, Howe/Hyslop 95% upper bound of
  (muT - muR)^2 - theta*sWR^2 <= 0, theta = (ln 1.25 / 0.25)^2; otherwise ABE.
`:ema_abel`: if CVwR > 30%, the 90% CI must lie within exp(+-0.760*sWR) with the
  sWR capped at CVwR 50% (69.84-143.19%); otherwise ABE.
Both also require the point estimate within 0.80-1.25.
"""
function scaled_be_replicate(t1, t2, r1, r2; method = :fda_rsabe, alpha = 0.05)
    n   = length(t1)
    df  = n - 2
    I   = (t1 .+ t2) ./ 2 .- (r1 .+ r2) ./ 2
    mu  = mean(I)
    se  = std(I) / sqrt(n)
    swr = sqrt(var(r1 .- r2) / 2)
    tq  = quantile(TDist(df), 1 - alpha)
    pe_ok = log(BE_LOWER) <= mu <= log(BE_UPPER)
    abe_ok() = (mu - tq * se >= log(BE_LOWER)) && (mu + tq * se <= log(BE_UPPER))
    if method == :fda_rsabe
        if swr < RSABE_SWR_MIN
            return (GMR = exp(mu), sWR = swr, scaled = false, bound = NaN,
                    BE = abe_ok() && pe_ok)
        end
        Em = mu^2
        Es = -RSABE_THETA * swr^2
        Cm = (abs(mu) + tq * se)^2
        Cs = Es * df / quantile(Chisq(df), 1 - alpha)
        bound = Em + Es + sqrt((Cm - Em)^2 + (Cs - Es)^2)
        return (GMR = exp(mu), sWR = swr, scaled = true, bound = bound,
                BE = (bound <= 0) && pe_ok)
    elseif method == :ema_abel
        cvwr = sqrt(exp(swr^2) - 1)
        cvwr <= 0.30 && return (GMR = exp(mu), sWR = swr, scaled = false, bound = NaN,
                                BE = abe_ok() && pe_ok)
        lim = ABEL_K * min(swr, ABEL_SWR_CAP)
        ok  = (mu - tq * se >= -lim) && (mu + tq * se <= lim)
        return (GMR = exp(mu), sWR = swr, scaled = true, bound = lim,
                BE = ok && pe_ok)
    else
        error("Unknown method $method (use :fda_rsabe or :ema_abel)")
    end
end

"""
    vbe_replicate_scaled(test, ref; n_subj, n_trial, cv_wr, cv_wt = cv_wr,
                         method = :fda_rsabe, seed) -> NamedTuple

Virtual full-replicate (4-period) trials: bootstrap `n_subj` subjects from the
PBPK pairs, give each two T and two R administrations with independent
within-subject noise (CV `cv_wt`, `cv_wr`), and apply `scaled_be_replicate`.
Returns the BE pass rate and the share of trials in which scaling applied.
"""
function vbe_replicate_scaled(test, ref; n_subj = nothing, n_trial = 2000,
                              cv_wr, cv_wt = cv_wr, method = :fda_rsabe, seed = 1234)
    rng  = MersenneTwister(seed)
    keep = @. !ismissing(test) & !ismissing(ref) &
              isfinite(test) & isfinite(ref) & (test > 0) & (ref > 0)
    lt0  = log.(Float64.(test[keep])); lr0 = log.(Float64.(ref[keep]))
    N    = length(lr0)
    nss  = n_subj === nothing ? N : n_subj
    swT  = sqrt(log(cv_wt^2 + 1)); swR = sqrt(log(cv_wr^2 + 1))
    npass = 0; nscaled = 0; gmrs = Vector{Float64}(undef, n_trial)
    for k in 1:n_trial
        idx = rand(rng, 1:N, nss)
        lt  = lt0[idx]; lr = lr0[idx]
        res = scaled_be_replicate(lt .+ swT .* randn(rng, nss), lt .+ swT .* randn(rng, nss),
                                  lr .+ swR .* randn(rng, nss), lr .+ swR .* randn(rng, nss);
                                  method)
        npass += res.BE; nscaled += res.scaled; gmrs[k] = res.GMR
    end
    (n_trial = n_trial, n_subj = nss, cv_wr = cv_wr, cv_wt = cv_wt, method = method,
     pass_rate = npass / n_trial, scaled_rate = nscaled / n_trial, GMR_med = median(gmrs))
end

end # module
