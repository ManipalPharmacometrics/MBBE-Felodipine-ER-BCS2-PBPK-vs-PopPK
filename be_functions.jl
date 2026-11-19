# =============================================================================
# be_functions.jl — the bioequivalence (BE) calculations used in steps 4–6.
#
#   be_test             one paired crossover comparison: GMR, 90% CI, pass/fail
#   pass_rate           share of simulated 2x2 BE trials that pass (one endpoint)
#   pass_rate_joint     share that pass on Cmax AND AUC together
#   pass_rate_scaled    share of simulated replicate trials passing FDA RSABE or EMA ABEL
#   values_at           per-subject endpoint anywhere along a design axis
#   first_crossing      point on an axis where the GMR first leaves 0.80–1.25
#
# Design: PK-Sim simulates the SAME 100 virtual subjects under every
# formulation, so test vs reference is paired (crossover). BE is tested on the
# within-subject log differences (two one-sided tests = 90% CI in 0.80–1.25).
# =============================================================================

using Distributions, Random, Statistics

"Paired average BE for one endpoint: GMR, 90% CI and pass (CI inside 0.80–1.25)."
function be_test(test, ref)
    d  = log.(test) .- log.(ref)
    n  = length(d)
    se = std(d) / sqrt(n)
    tq = quantile(TDist(n - 1), 0.95)
    lo, hi = exp(mean(d) - tq * se), exp(mean(d) + tq * se)
    (n = n, GMR = exp(mean(d)), CI_lower = lo, CI_upper = hi, pass = lo >= BE_LOWER && hi <= BE_UPPER)
end

"""
    pass_rate(test, ref; n_subj, n_trial, cv)

Simulate `n_trial` 2x2 BE trials: draw `n_subj` subjects with replacement,
add within-subject noise (CV `cv`) to each test and reference value, run
`be_test`. Returns the fraction that pass (power, or type I error at a boundary).
"""
function pass_rate(test, ref; n_subj, n_trial, cv)
    rng = MersenneTwister(SEED_TRIALS)
    sw  = sqrt(log(cv^2 + 1))
    N   = length(ref)
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, n_subj)
        tv  = test[idx] .* exp.(randn(rng, n_subj) .* sw)
        rv  = ref[idx]  .* exp.(randn(rng, n_subj) .* sw)
        npass += be_test(tv, rv).pass
    end
    npass / n_trial
end

"As `pass_rate`, but a trial passes only if Cmax AND AUC both pass (same subjects)."
function pass_rate_joint(testC, refC, testA, refA; n_subj, n_trial, cvC, cvA)
    rng = MersenneTwister(SEED_TRIALS)
    swC = sqrt(log(cvC^2 + 1)); swA = sqrt(log(cvA^2 + 1))
    N   = length(refC)
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, n_subj)
        tc = testC[idx] .* exp.(randn(rng, n_subj) .* swC)
        rc = refC[idx]  .* exp.(randn(rng, n_subj) .* swC)
        ta = testA[idx] .* exp.(randn(rng, n_subj) .* swA)
        ra = refA[idx]  .* exp.(randn(rng, n_subj) .* swA)
        npass += be_test(tc, rc).pass && be_test(ta, ra).pass
    end
    npass / n_trial
end

"""
    scaled_trial(t1, t2, r1, r2, method)

One full-replicate trial (each subject gets test twice and reference twice),
on log values. sWR = reference within-subject SD.
  :fda_rsabe — if sWR >= 0.294: upper 95% bound of (muT-muR)^2 - θ·sWR^2 <= 0,
               θ = (ln 1.25 / 0.25)^2 (Howe's method); else plain ABE.
  :ema_abel  — if CVwR > 30%: 90% CI within exp(±0.760·sWR), sWR capped at
               CVwR 50%; else plain ABE.
Both also need the point estimate inside 0.80–1.25.
"""
function scaled_trial(t1, t2, r1, r2, method)
    n  = length(t1); df = n - 2
    I  = (t1 .+ t2) ./ 2 .- (r1 .+ r2) ./ 2
    mu = mean(I); se = std(I) / sqrt(n)
    swr = sqrt(var(r1 .- r2) / 2)
    tq  = quantile(TDist(df), 0.95)
    pe_ok  = log(BE_LOWER) <= mu <= log(BE_UPPER)
    abe_ok = mu - tq * se >= log(BE_LOWER) && mu + tq * se <= log(BE_UPPER)
    if method == :fda_rsabe
        swr < 0.294 && return abe_ok && pe_ok
        θ  = (log(1.25) / 0.25)^2
        Em, Es = mu^2, -θ * swr^2
        Cm = (abs(mu) + tq * se)^2
        Cs = Es * df / quantile(Chisq(df), 0.95)
        return Em + Es + sqrt((Cm - Em)^2 + (Cs - Es)^2) <= 0 && pe_ok
    else  # :ema_abel
        sqrt(exp(swr^2) - 1) <= 0.30 && return abe_ok && pe_ok
        lim = 0.760 * min(swr, sqrt(log(0.50^2 + 1)))
        return mu - tq * se >= -lim && mu + tq * se <= lim && pe_ok
    end
end

"Share of simulated replicate trials (Cmax) that pass RSABE or ABEL."
function pass_rate_scaled(test, ref; n_subj, n_trial, method)
    rng = MersenneTwister(SEED_TRIALS)
    sw  = sqrt(log(CV_WR^2 + 1))                 # same CV for test and reference
    lt0, lr0 = log.(test), log.(ref)
    N = length(lr0)
    npass = 0
    for _ in 1:n_trial
        idx = rand(rng, 1:N, n_subj)
        lt, lr = lt0[idx], lr0[idx]
        npass += scaled_trial(lt .+ sw .* randn(rng, n_subj), lt .+ sw .* randn(rng, n_subj),
                              lr .+ sw .* randn(rng, n_subj), lr .+ sw .* randn(rng, n_subj), method)
    end
    npass / n_trial
end

# ---- moving along one axis of the design ------------------------------------

"Per-subject endpoint of one formulation, ordered by IndividualId."
subject_values(pk, t50, shape, endpoint) =
    sort(pk[(pk.T50 .== t50) .& (pk.shape .== shape), :], :IndividualId)[!, endpoint]

"Formulations on an axis (the values of T50 or shape), sorted."
axis_points(pk, ax) = sort(unique(pk[pk[!, ax.fixed] .== ax.fixed_value, ax.col]))

"Endpoint values of the formulation at position `x` on axis `ax` (a simulated one)."
function at_point(pk, ax, x, endpoint)
    ax.name == "T50" ? subject_values(pk, x, ax.fixed_value, endpoint) :
                       subject_values(pk, ax.fixed_value, x, endpoint)
end

"""
    values_at(pk, ax, endpoint, x)

Per-subject `endpoint` at any position `x` on axis `ax`, interpolating each
subject's log(value) between the two neighbouring formulations (log(T50) on the
T50 axis, shape on the shape axis).
"""
function values_at(pk, ax, endpoint, x)
    A = axis_points(pk, ax)
    x = clamp(x, A[1], A[end])
    x in A && return at_point(pk, ax, x, endpoint)
    k = searchsortedlast(A, x)
    w = ax.log_x ? (log(x) - log(A[k])) / (log(A[k+1]) - log(A[k])) : (x - A[k]) / (A[k+1] - A[k])
    exp.((1 - w) .* log.(at_point(pk, ax, A[k], endpoint)) .+ w .* log.(at_point(pk, ax, A[k+1], endpoint)))
end

"161-point grid along an axis (log-spaced for T50)."
grid_of(pk, ax) = (A = axis_points(pk, ax);
                   ax.log_x ? exp.(range(log(A[1]), log(A[end]), length = GRID_POINTS)) :
                              collect(range(A[1], A[end], length = GRID_POINTS)))

"True (noise-free) GMR at position `x`."
gmr_at(pk, ax, endpoint, x) = be_test(values_at(pk, ax, endpoint, x), values_at(pk, ax, endpoint, ax.ref)).GMR

"Pass rate of 2x2 trials at position `x` (endpoint Cmax, AUC_tEnd or Joint)."
function pass_rate_at(pk, ax, endpoint, x, n_subj, n_trial)
    if endpoint == "Joint"
        pass_rate_joint(values_at(pk, ax, "Cmax", x), values_at(pk, ax, "Cmax", ax.ref),
                        values_at(pk, ax, "AUC_tEnd", x), values_at(pk, ax, "AUC_tEnd", ax.ref);
                        n_subj, n_trial, cvC = WSV_CV["Cmax"], cvA = WSV_CV["AUC_tEnd"])
    else
        pass_rate(values_at(pk, ax, endpoint, x), values_at(pk, ax, endpoint, ax.ref);
                  n_subj, n_trial, cv = WSV_CV[endpoint])
    end
end

"""
    first_crossing(grid, gmr, ref, side)

Walk from the reference down (`:low`) or up (`:high`) the grid and return the
first position where the GMR crosses 0.80 or 1.25 (linear interpolation), and
which limit. `missing` if it never does inside the simulated range.
"""
function first_crossing(grid, gmr, ref, side)
    i = argmin(abs.(grid .- ref))
    step = side === :low ? -1 : 1
    while 1 <= i + step <= length(grid)
        j = i + step
        for (name, limit) in (("0.80", BE_LOWER), ("1.25", BE_UPPER))
            a, b = gmr[i], gmr[j]
            if (a - limit) * (b - limit) <= 0 && a != b
                w = (limit - a) / (b - a)
                return (x = grid[i] + w * (grid[j] - grid[i]), limit = name)
            end
        end
        i = j
    end
    (x = missing, limit = "not reached")
end
