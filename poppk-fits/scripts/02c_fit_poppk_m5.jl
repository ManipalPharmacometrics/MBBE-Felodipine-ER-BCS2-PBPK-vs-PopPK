# =============================================================================
# 02c_fit_poppk_m5.jl — as 02b_fit_poppk_v2.jl, with the final model M5 (dual
# Weibull @delay + release lag -> ka -> 2-cmt; src/poppk_models_felo_structure.jl),
# 2026-10-24. M3 (output_v2/) is kept as the fallback.
# Default output: output_m5/ (be_dataset.csv copied from output/ if absent).
#
# Starts as 06_structure_fits.jl, keep higher LL (default one start, see STARTS). Fits already done in
# output_modelcomp/M5/<arm>/<start>/fit.jls (the 5 diagnostic arms) are reused.
# params.csv also holds summaries of the mixture release so felo_05 runs as
# for M3: tvt50 = time of 50% released (lag excluded), tvgamma = shape of the
# single Weibull through the mixture's 25% and 75% release times.
# Resumable: one folder per arm under output_m5/arms/, skipped if done.
# Optional env: POPPK_ARMS, POPPK_SHARD as 02b; POPPK_STARTS, POPPK_REUSE (below).
# Run: julia -J/opt/julia/lib/julia/pumas.so --project=<repo>/julia-env scripts/02c_fit_poppk_m5.jl
# =============================================================================

using Pumas, CSV, DataFrames, Statistics, Serialization

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output_m5")); mkpath(OUT)
const ARMS = joinpath(OUT, "arms"); mkpath(ARMS)
const MC   = joinpath(ROOT, "output_modelcomp", "M5")
module Smod; using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_models_felo_structure.jl")); end   # module name as 06, so its fit.jls deserialize
const poppk_m5 = Smod.poppk_m5
for f in ("be_dataset.csv", "pbpk_pk_parameters_individual.csv")
    isfile(joinpath(OUT, f)) || cp(joinpath(ROOT, "output", f), joinpath(OUT, f))
end

const DENSE = collect(0.0:0.05:96.0)
data = CSV.read(joinpath(OUT, "be_dataset.csv"), DataFrame)
arms = sort(unique(data[:, [:T50, :shape]]), [:shape, :T50])
key(t, s) = "$(t)_$(s)"
if haskey(ENV, "POPPK_ARMS")
    want = Set(split(ENV["POPPK_ARMS"], ","))
    arms = filter(r -> key(r.T50, r.shape) in want, arms)
end
if haskey(ENV, "POPPK_SHARD")
    k, m = parse.(Int, split(ENV["POPPK_SHARD"], "/"))
    arms = arms[k:m:nrow(arms), :]
end
# arms with reusable diagnostic fits first (quick check of the output code)
arms = arms[sortperm([!isdir(joinpath(MC, key(r.T50, r.shape))) for r in eachrow(arms)]), :]

trap(t, c) = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in 2:length(t))

# starts as 06_structure_fits.jl
starts(t, s) = [
    "release_limited"    => (tvt50a = 0.5t, tvga = s, tvt50b = 3t, tvgb = s, tvff = 0.6, tvka = 2.0,
                             tvcl = 300.0, tvvc = 200.0, tvq = 200.0, tvvp = 1500.0, tvlag = 0.3113),
    "absorption_limited" => (tvt50a = 0.3t, tvga = 1.0, tvt50b = 2t, tvgb = 1.0, tvff = 0.3, tvka = 0.5,
                             tvcl = 300.0, tvvc = 500.0, tvq = 100.0, tvvp = 1000.0, tvlag = 0.3113)]

# Start(s) to run, env POPPK_STARTS (default absorption-limited only): on the 7
# arms fitted with both starts (5 diagnostic + 4.07/4.78 x 0.6) it was never
# worse than release-limited and the higher LL on 3, so the second start (half
# the run time, ~30 min per arm single-threaded) is skipped.
const STARTS = split(get(ENV, "POPPK_STARTS", "absorption_limited"), ",")
# Reuse the diagnostic fits only for the default dataset: POPPK_REUSE = "0" to refit
const REUSE = get(ENV, "POPPK_REUSE", "1") == "1"

# mixture release fraction F(t) (lag excluded) and the time it reaches p (bisection)
relF(c, t) = (wb(t50, g) = 1 - 2.0^(-(t / t50)^g); c.tvff * wb(c.tvt50a, c.tvga) + (1 - c.tvff) * wb(c.tvt50b, c.tvgb))
function tq(c, p)
    lo, hi = 0.0, 1.0
    while relF(c, hi) < p; hi *= 2; end
    for _ in 1:60; mid = (lo + hi) / 2; relF(c, mid) < p ? (lo = mid) : (hi = mid); end
    (lo + hi) / 2
end

for a in eachrow(arms)
    dir = joinpath(ARMS, key(a.T50, a.shape))
    isfile(joinpath(dir, "nca.csv")) && (println("skip ", key(a.T50, a.shape)); continue)
    mkpath(dir)
    arm = filter(r -> r.T50 == a.T50 && r.shape == a.shape, data)
    arm.cmt .= 1                                # Gut, bioav 0 (input via @delay)
    pop = read_pumas(arm; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                     observations = [:dv], covariates = [:T50, :shape, :IndividualId])
    @info "Fitting T50 = $(a.T50) h, shape = $(a.shape) (n = $(length(pop)))"
    t0 = time()
    fits = []
    for (lab, p) in filter(x -> x[1] in STARTS, starts(a.T50, a.shape))
        jls = joinpath(MC, key(a.T50, a.shape), lab, "fit.jls")
        REUSE || (jls = "")                     # other datasets (e.g. noise replicates) must refit
        f = try
            isfile(jls) ? deserialize(jls) :
            fit(poppk_m5, pop, merge(init_params(poppk_m5), p), FOCE();
                constantcoef = (:tvlag,), optim_options = (show_trace = false,))
        catch err
            @warn "fit failed" lab err; nothing
        end
        f === nothing || push!(fits, (lab = lab, f = f, ll = loglikelihood(f)))
    end
    isempty(fits) && (@warn "no fit for arm" a.T50 a.shape; continue)
    best = fits[argmax([x.ll for x in fits])]
    fpm = best.f
    el = round(time() - t0, digits = 1)
    c = coef(fpm)
    ofv = try; deviance(fpm); catch; missing; end
    t25, t50, t75 = tq(c, 0.25), tq(c, 0.5), tq(c, 0.75)
    CSV.write(joinpath(dir, "params.csv"),
        DataFrame(T50 = a.T50, shape = a.shape, tvt50 = t50,
                  tvgamma = log(log(4) / log(4 / 3)) / log(t75 / t25),
                  tvt50a = c.tvt50a, tvga = c.tvga, tvt50b = c.tvt50b, tvgb = c.tvgb, tvff = c.tvff,
                  tvka = c.tvka, tvcl = c.tvcl, tvvc = c.tvvc, tvq = c.tvq, tvvp = c.tvvp, tvlag = c.tvlag,
                  omega_t50 = c.Ω[1, 1], omega_cl = c.Ω[2, 2], omega_vc = c.Ω[3, 3],
                  sigma_p = c.σ_p, OFV = ofv, LL = best.ll, basin = best.lab,
                  LL_other = length(fits) == 2 ? minimum(x.ll for x in fits) : missing,
                  fit_seconds = el))
    try
        ct = DataFrame(coeftable(infer(fpm))); insertcols!(ct, 1, :T50 => a.T50, :shape => a.shape)
        CSV.write(joinpath(dir, "se.csv"), ct)
    catch err
        @warn "infer() failed" a.T50 a.shape err
    end
    ηs = ηshrinkage(fpm).η
    CSV.write(joinpath(dir, "shrinkage.csv"),
        DataFrame(T50 = a.T50, shape = a.shape, eta_t50 = ηs[1], eta_cl = ηs[2], eta_vc = ηs[3]))
    try; serialize(joinpath(dir, "fit.jls"), fpm); catch err; @warn "serialize failed" err; end
    preds = predict(fpm; obstimes = DENSE)
    nca = DataFrame(); thin = DataFrame()
    for (subj, pr) in zip(pop, preds)
        id = Int(subj.covariates(0.0).IndividualId)
        tt = pr.time; ip = pr.ipred.dv
        k = argmax(ip)
        push!(nca, (T50 = a.T50, shape = a.shape, IndividualId = id, Cmax = ip[k],
                    Tmax = tt[k], AUC_tEnd = trap(tt, ip)))
        for j in 1:10:length(tt)
            push!(thin, (T50 = a.T50, shape = a.shape, IndividualId = id, time = tt[j], ipred = ip[j]))
        end
    end
    CSV.write(joinpath(dir, "ipred_thin.csv"), thin)
    CSV.write(joinpath(dir, "nca.csv"), nca)          # written last = arm done
    println("done $(key(a.T50, a.shape)) [$(best.lab)] in $(el) s  T50eff=$(round(t50, digits=2))  ff=$(round(c.tvff, digits=2))  CL=$(round(c.tvcl, digits=1))  sigma=$(round(c.σ_p, digits=3))"); flush(stdout)
end
