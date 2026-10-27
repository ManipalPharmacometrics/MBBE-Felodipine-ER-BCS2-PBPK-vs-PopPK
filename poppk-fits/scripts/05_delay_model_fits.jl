# =============================================================================
# 05_delay_model_fits.jl — refit M1, M2, M3 in the @delay form
# (src/poppk_models_felo_delay.jl) on the same 5 diagnostic arms and starts as
# 04, saving every fit, so the fitted results can be compared with the
# hand-coded Weibull ODE versions (2026-10-24).
#
# Variants (release lag held constant, as in 04):
#   M1d, M2d, M3d     lag 0.3113 h (same as M1–M3)
#   M1d0, M2d0, M3d0  lag 0 (the plain @delay form has no lag)
# Same outputs per fit as 04, in output_modelcomp/<model>/<arm>/<start>/.
# All records cmt 1 (dose compartment with bioav 0). M4 is not refitted: it
# was rejected, and its 200 h period-2 dose would need the same truncation.
# Resumable: a fit with fitstats.csv is skipped.
# Env: MC_MODELS = "M1d,M2d,M3d,M1d0,M2d0,M3d0" (default all), MC_ARMS = "4.78_0.5,..."
# Run: julia -J/opt/julia/lib/julia/pumas.so --project=<repo>/julia-env scripts/05_delay_model_fits.jl
# =============================================================================

using Pumas, CSV, DataFrames, Statistics, LinearAlgebra, Serialization

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = joinpath(ROOT, "output_modelcomp"); mkpath(OUT)
module M1mod; using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_model_felo_1cmt.jl")); end
module M3mod; using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_model_felo_v2.jl")); end
module Dmod;  using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_models_felo_delay.jl")); end

const REF = (4.78, 0.764)
const ARMS_ALL = [(4.78, 0.764), (4.78, 0.5), (4.78, 1.5), (2.39, 0.764), (9.56, 0.764)]
akey(t, s) = "$(t)_$(s)"

# starts: M1 and M3/M4 as in production; M2 = the two 10-16 basins
const M2_INITS = [
    "release_limited"     => (tvalpha = 7.0, tvgamma = 0.9, tvcl = 300.0, tvvc = 250.0, tvq = 200.0, tvvp = 1500.0),
    "disposition_limited" => (tvalpha = 1.0, tvgamma = 1.3, tvcl = 300.0, tvvc = 2000.0, tvq = 300.0, tvvp = 3000.0)]
# estimated parameters (tvlag fixed in every model): θ + Ω diagonal + σ
base(m) = m[1:3]                                        # "M3d0" -> "M3d"
const NPAR = Dict("M1d" => 4 + 3 + 1, "M2d" => 6 + 3 + 1, "M3d" => 7 + 3 + 1)
const MODEL = Dict("M1d" => Dmod.poppk_m1d, "M2d" => Dmod.poppk_m2d, "M3d" => Dmod.poppk_m3d)
lagof(m) = endswith(m, "d0") ? 0.0 : 0.3113
TP2 = M3mod.TP2

data = CSV.read(joinpath(ROOT, "output", "be_dataset.csv"), DataFrame)
pb   = CSV.read(joinpath(ROOT, "output", "pbpk_pk_parameters_individual.csv"), DataFrame)
armdata(t, s) = filter(r -> r.T50 == t && r.shape == s, data)
pbarm(t, s)   = sort(filter(r -> r.T50 == t && r.shape == s, pb), :IndividualId)
trap(t, c) = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in 2:length(t))

function population(m, t, s)
    d = armdata(t, s); d.cmt .= 1                       # Central (M1d/M2d) or Gut (M3d)
    read_pumas(d; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
               observations = [:dv], covariates = [:IndividualId])
end

function starts(m, t, s)
    lag = (tvlag = lagof(m),)
    base(m) == "M1d" && return [l => merge(p, lag) for (l, p) in M1mod.INITS_FELO]
    base(m) == "M2d" && return [l => merge(p, lag) for (l, p) in M2_INITS]
    [l => merge(p, (tvt50 = t, tvgamma = s), lag) for (l, p) in M3mod.INITS_V2]
end

function coefrows(c)
    rows = DataFrame(param = String[], value = Float64[])
    for (k, v) in pairs(c)
        if v isa AbstractMatrix
            for i in 1:size(v, 1); push!(rows, ("$(k)[$i,$i]", v[i, i])); end
        else
            push!(rows, (string(k), v))
        end
    end
    rows
end

function endpoints(m, pop, f, t, s)
    grid = m == "M4" ? vcat(collect(0.0:0.05:96.0), collect(TP2:0.05:TP2 + 96.0)) : collect(0.0:0.05:96.0)
    pr = predict(f; obstimes = grid)
    ep = DataFrame(); thin = DataFrame()
    for (subj, p) in zip(pop, pr)
        id = Int(subj.covariates(0.0).IndividualId); tt = p.time; ip = p.ipred.dv
        periods = m == "M4" ? [("R", REF..., tt .<= 96.0, 0.0), ("T", t, s, tt .>= TP2, TP2)] :
                              [("T", t, s, trues(length(tt)), 0.0)]
        for (per, T50, sh, sel, off) in periods
            ts = tt[sel] .- off; cs = ip[sel]; k = argmax(cs)
            push!(ep, (period = per, T50 = T50, shape = sh, IndividualId = id,
                       Cmax = cs[k], Tmax = ts[k], AUC_tEnd = trap(ts, cs)))
            for j in 1:10:length(ts)
                push!(thin, (period = per, T50 = T50, shape = sh, IndividualId = id, time = ts[j], ipred = cs[j]))
            end
        end
    end
    ep = leftjoin(ep, select(pb, :T50, :shape, :IndividualId, :Cmax => :Cmax_PBPK, :AUC_tEnd => :AUC_PBPK),
                  on = [:T50, :shape, :IndividualId])
    sort!(ep, [:period, :IndividualId]), thin
end

want_m = split(get(ENV, "MC_MODELS", "M1d,M2d,M3d,M1d0,M2d0,M3d0"), ",")
want_a = haskey(ENV, "MC_ARMS") ? Set(split(ENV["MC_ARMS"], ",")) : nothing
for m in want_m, (t, s) in ARMS_ALL
    want_a === nothing || akey(t, s) in want_a || continue
    pop = population(m, t, s); model = MODEL[base(m)]
    nobs = sum(count(!ismissing, sb.observations.dv) for sb in pop)
    for (lab, p) in starts(m, t, s)
        dir = joinpath(OUT, m, akey(t, s), lab)
        isfile(joinpath(dir, "fitstats.csv")) && (println("skip $m $(akey(t, s)) $lab"); continue)
        mkpath(dir); t0 = time()
        f = try
            fit(model, pop, merge(init_params(model), p), FOCE();
                constantcoef = (:tvlag,), optim_options = (show_trace = false,))
        catch err
            @warn "fit failed" m t s lab err; err
        end
        secs = round(time() - t0, digits = 1)
        if f isa Exception
            CSV.write(joinpath(dir, "fitstats.csv"),
                DataFrame(model = m, arm = akey(t, s), start = lab, status = "failed: " * first(sprint(showerror, f), 200),
                          LL = missing, m2LL = missing, n_par = NPAR[base(m)], n_obs = nobs, AIC = missing, BIC = missing,
                          sigma = missing, secs = secs))
            continue
        end
        ll = loglikelihood(f); k = NPAR[base(m)]
        status = try; Pumas.Optim.converged(f.optim) ? "converged" : "not converged"; catch; "unknown"; end
        CSV.write(joinpath(dir, "coef.csv"), coefrows(coef(f)))
        try
            CSV.write(joinpath(dir, "se.csv"), DataFrame(coeftable(infer(f))))
        catch err
            @warn "infer failed" m t s lab
        end
        try
            CSV.write(joinpath(dir, "shrinkage.csv"), DataFrame(eta = 1:length(ηshrinkage(f).η), shrinkage = collect(ηshrinkage(f).η)))
        catch; end
        ep, thin = endpoints(m, pop, f, t, s)
        CSV.write(joinpath(dir, "endpoints.csv"), ep); CSV.write(joinpath(dir, "ipred_thin.csv"), thin)
        try; serialize(joinpath(dir, "fit.jls"), f); catch err; @warn "serialize failed" err; end
        CSV.write(joinpath(dir, "fitstats.csv"),           # written last = fit done
            DataFrame(model = m, arm = akey(t, s), start = lab, status = status, LL = ll, m2LL = -2ll,
                      n_par = k, n_obs = nobs, AIC = -2ll + 2k, BIC = -2ll + k * log(nobs),
                      sigma = coef(f).σ_p, secs = secs))
        println("done $m $(akey(t, s)) $lab  LL=$(round(ll, digits=1))  [$status, $(secs) s]"); flush(stdout)
    end
end
