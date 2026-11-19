# =============================================================================
# STEP 6b — cross-check: do our fast BE tests give the same pass/fail as the
# standard Pumas test, `pumas_be`?
#
# In : results/step1_pbpk_endpoints.csv, results/step3_poppk_endpoints.csv,
#      results/step6_boundaries.csv
# Out: results/step6b_check_vs_pumas_be.csv
#
# Why: steps 5–6 simulate hundreds of thousands of BE trials, so they use the
# fast tests in be_functions.jl: `be_test` (2x2 ABE) and `scaled_trial`
# (FDA RSABE, EMA ABEL). Step 4 uses `pumas_be`. This step shows they reach
# the same decisions.
#
# Where: at each BE boundary found in step 6, where the true GMR is 0.80 or
# 1.25, so trials sit as close to the pass/fail line as they get.
#   ABE   : Cmax and AUC, both axes, both sides, n = 24 and 54 (2x2, RT/TR)
#   RSABE, ABEL : Cmax boundary, n = 24 (full replicate, RTRT/TRTR), one
#                 pumas_be call per trial with the FDA_HighlyVariable criteria.
#     RSABE = pumas_be's own pass/fail (`assess_be`).
#     ABEL  = the EMA limits applied to pumas_be's 90% CI and CVwR:
#             CVwR > 30%: CI within exp(±0.760·sWR), sWR capped at CVwR 50%;
#             otherwise 0.80–1.25; and the GMR within 0.80–1.25.
#     Not `EMA_HighlyVariable`: in this Bioequivalence version it runs the FDA
#     progesterone RSABE with θ = 0.76, without the EMA widened limits or the
#     50% cap (pumas_be prints a warning saying so), so it is not EMA ABEL.
#
# How: the first N_TRIAL_CHECK trials of `pass_rate` / `pass_rate_scaled`
# (same seed, same draws), each one run through BOTH tests.
# =============================================================================

using CSV, DataFrames, Random, Statistics, Bioequivalence
include("settings.jl")
include("be_functions.jl")

pbpk  = CSV.read(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), DataFrame)
poppk = CSV.read(joinpath(RESULTS, "step3_poppk_endpoints.csv"), DataFrame)
bnd   = CSV.read(joinpath(RESULTS, "step6_boundaries.csv"), DataFrame)
SOURCES = Dict("PBPK" => pbpk, "PopPK" => poppk)
out = DataFrame()

# ---- ABE, 2x2 crossover ---------------------------------------------------------
for b in eachrow(bnd), n in N_SUBJ
    (b.endpoint in ENDPOINTS && !ismissing(b.boundary)) || continue
    ax = only(a for a in AXES if a.name == b.axis)
    pk = SOURCES[b.source]
    test, ref = values_at(pk, ax, b.endpoint, b.boundary), values_at(pk, ax, b.endpoint, ax.ref)
    sw  = sqrt(log(WSV_CV[b.endpoint]^2 + 1))
    rng = MersenneTwister(SEED_TRIALS)
    seq = [isodd(i) ? "RT" : "TR" for i in 1:n]
    crossover = DataFrame(id = repeat(1:n, inner = 2), sequence = repeat(seq, inner = 2),
                          period = repeat([1, 2], n))
    pass_own, pass_pumas = falses(N_TRIAL_CHECK), falses(N_TRIAL_CHECK)
    for k in 1:N_TRIAL_CHECK
        idx = rand(rng, 1:length(ref), n)                           # as in pass_rate
        tv  = test[idx] .* exp.(randn(rng, n) .* sw)
        rv  = ref[idx]  .* exp.(randn(rng, n) .* sw)
        pass_own[k] = be_test(tv, rv).pass
        crossover.value = [s[p] == 'T' ? tv[i] : rv[i]
                           for (i, s, p) in zip(crossover.id, crossover.sequence, crossover.period)]
        r = pumas_be(crossover; endpoint = :value, sigdigits = 16).result[1, :]
        pass_pumas[k] = r.LB >= BE_LOWER && r.UB <= BE_UPPER
    end
    push!(out, (method = "ABE", axis = b.axis, source = b.source, endpoint = b.endpoint, side = b.side,
                boundary = b.boundary, n_subj = n, n_trial = N_TRIAL_CHECK,
                pass_rate_ours = mean(pass_own), pass_rate_pumas_be = mean(pass_pumas),
                same_decision = mean(pass_own .== pass_pumas)))
end

# ---- RSABE / ABEL, full replicate (Cmax) -------------------------------------------
sw  = sqrt(log(CV_WR^2 + 1))
seq = [isodd(i) ? "RTRT" : "TRTR" for i in 1:N_SCALED]
replicate = DataFrame(id = repeat(1:N_SCALED, inner = 4), sequence = repeat(seq, inner = 4),
                      period = repeat(1:4, N_SCALED))
for b in eachrow(bnd[(bnd.endpoint .== "Cmax") .& .!ismissing.(bnd.boundary), :])
    ax = only(a for a in AXES if a.name == b.axis)
    pk = SOURCES[b.source]
    lt0 = log.(values_at(pk, ax, "Cmax", b.boundary)); lr0 = log.(values_at(pk, ax, "Cmax", ax.ref))
    rng = MersenneTwister(SEED_TRIALS)
    own   = Dict(m => falses(N_TRIAL_CHECK) for m in (:fda_rsabe, :ema_abel))
    pumas = Dict(m => falses(N_TRIAL_CHECK) for m in (:fda_rsabe, :ema_abel))
    for k in 1:N_TRIAL_CHECK
        idx = rand(rng, 1:length(lr0), N_SCALED)                    # as in pass_rate_scaled
        lt, lr = lt0[idx], lr0[idx]
        t1 = lt .+ sw .* randn(rng, N_SCALED); t2 = lt .+ sw .* randn(rng, N_SCALED)
        r1 = lr .+ sw .* randn(rng, N_SCALED); r2 = lr .+ sw .* randn(rng, N_SCALED)
        own[:fda_rsabe][k] = scaled_trial(t1, t2, r1, r2, :fda_rsabe)
        own[:ema_abel][k]  = scaled_trial(t1, t2, r1, r2, :ema_abel)
        # 1st T / 1st R go to the earlier period, 2nd T / 2nd R to the later one
        replicate.value = [s[p] == 'T' ? (p <= 2 ? t1[i] : t2[i]) : (p <= 2 ? r1[i] : r2[i])
                           for (i, s, p) in zip(replicate.id, replicate.sequence, replicate.period)]
        o = pumas_be(replicate, FDA_HighlyVariable; endpoint = :value, logtransformed = true, sigdigits = 16)
        pumas[:fda_rsabe][k] = all(assess_be(FDA_HighlyVariable, o).assessments)
        r   = o.result[1, :]
        lim = r.CVᵣ > 0.30 ? exp(0.760 * sqrt(log(min(r.CVᵣ, 0.50)^2 + 1))) : BE_UPPER
        pumas[:ema_abel][k] = r.LB >= 1 / lim && r.UB <= lim && BE_LOWER <= r.GMR <= BE_UPPER
    end
    for (m, label) in (:fda_rsabe => "RSABE (FDA)", :ema_abel => "ABEL (EMA)")
        push!(out, (method = label, axis = b.axis, source = b.source, endpoint = "Cmax", side = b.side,
                    boundary = b.boundary, n_subj = N_SCALED, n_trial = N_TRIAL_CHECK,
                    pass_rate_ours = mean(own[m]), pass_rate_pumas_be = mean(pumas[m]),
                    same_decision = mean(own[m] .== pumas[m])))
    end
end

CSV.write(joinpath(RESULTS, "step6b_check_vs_pumas_be.csv"), out)
println("Pass rate at the BE boundary, our tests vs pumas_be (same trials):")
show(out; allrows = true, allcols = true); println()
for m in unique(out.method)
    d = out[out.method .== m, :]
    println(rpad(m, 12), " lowest agreement ", round(100 * minimum(d.same_decision); digits = 2),
            "%, largest pass-rate difference ",
            round(100 * maximum(abs.(d.pass_rate_ours .- d.pass_rate_pumas_be)); digits = 2), " points")
end
