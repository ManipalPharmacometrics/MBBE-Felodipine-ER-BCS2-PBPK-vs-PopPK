# =============================================================================
# 04_model_comparison_fits.jl — fit every candidate felodipine PopPK structure
# on the same diagnostic arms, from every start, and save ALL fits (not only
# the best) so the models can be compared side by side (2026-10-23).
#
# Models (equations: Felodipine-BCSII/poppk/model_comparison/MODELS.md)
#   M1  1-cmt + Weibull release                   src/poppk_model_felo_1cmt.jl (v1 baseline)
#   M2  2-cmt + Weibull release                   src/poppk_model_felo.jl      (10-16, rejected)
#   M3  Weibull release -> ka -> 2-cmt, per arm   src/poppk_model_felo_v2.jl   (v2)
#   M4  as M3, ref + test paired, shared disposition  src/poppk_model_felo_v2.jl (paired)
# Arms: reference 4.78/0.764 and the corners 4.78/0.5, 4.78/1.5, 2.39/0.764,
# 9.56/0.764 (M4: the four non-reference arms, each paired with the reference).
#
# Per fit, output_modelcomp/<model>/<arm>/<start>/:
#   fitstats.csv  LL, -2LL, n_par, n_obs, AIC, BIC, sigma, runtime, status
#   coef.csv      every parameter (Ω diagonal expanded), estimated or fixed
#   se.csv        coeftable(infer(fit)) when infer succeeds
#   shrinkage.csv η-shrinkage
#   endpoints.csv per subject IPRED Cmax, AUC_0-96 vs PK-Sim (per period for M4)
#   ipred_thin.csv 0.5-h IPRED grid
#   fit.jls       serialized fit (reload in this env with Serialization)
# Resumable: a fit with fitstats.csv is skipped.
# Env: MC_MODELS = "M1,M2,M3,M4" (default all), MC_ARMS = "4.78_0.5,..."
# Run: julia -J/opt/julia/lib/julia/pumas.so --project=<repo>/julia-env scripts/04_model_comparison_fits.jl
# =============================================================================

using Pumas, CSV, DataFrames, Statistics, LinearAlgebra, Serialization

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = joinpath(ROOT, "output_modelcomp"); mkpath(OUT)
module M1mod; using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_model_felo_1cmt.jl")); end
module M2mod; using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_model_felo.jl")); end
module M3mod; using Pumas; include(joinpath(@__DIR__, "..", "src", "poppk_model_felo_v2.jl")); end

const REF = (4.78, 0.764)
const ARMS_ALL = [(4.78, 0.764), (4.78, 0.5), (4.78, 1.5), (2.39, 0.764), (9.56, 0.764)]
akey(t, s) = "$(t)_$(s)"

# starts: M1 and M3/M4 as in production; M2 = the two 10-16 basins
const M2_INITS = [
    "release_limited"     => (tvalpha = 7.0, tvgamma = 0.9, tvcl = 300.0, tvvc = 250.0, tvq = 200.0, tvvp = 1500.0),
    "disposition_limited" => (tvalpha = 1.0, tvgamma = 1.3, tvcl = 300.0, tvvc = 2000.0, tvq = 300.0, tvvp = 3000.0)]
# estimated parameters (tvlag fixed in every model): θ + Ω diagonal + σ
const NPAR = Dict("M1" => 4 + 3 + 1, "M2" => 6 + 3 + 1, "M3" => 7 + 3 + 1, "M4" => 9 + 3 + 1)
const MODEL = Dict("M1" => M1mod.poppk_model_felo, "M2" => M2mod.poppk_model_felo,
                   "M3" => M3mod.poppk_felo_wka2, "M4" => M3mod.poppk_felo_pair)
TP2 = M3mod.TP2

data = CSV.read(joinpath(ROOT, "output", "be_dataset.csv"), DataFrame)
pb   = CSV.read(joinpath(ROOT, "output", "pbpk_pk_parameters_individual.csv"), DataFrame)
armdata(t, s) = filter(r -> r.T50 == t && r.shape == s, data)
pbarm(t, s)   = sort(filter(r -> r.T50 == t && r.shape == s, pb), :IndividualId)
trap(t, c) = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in 2:length(t))

function population(m, t, s)
    if m == "M4"
        r = armdata(REF...); x = armdata(t, s)
        r.cmt = ifelse.(r.evid .== 1, 1, 4)            # Depot1 / Central
        x.cmt = ifelse.(x.evid .== 1, 2, 4); x.time .+= TP2
        d = vcat(r, x); d.ID = d.IndividualId; sort!(d, [:ID, :time])   # id column can't double as a covariate
        return read_pumas(d; id = :ID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                          observations = [:dv], covariates = [:IndividualId])
    end
    d = armdata(t, s)
    m == "M3" && (d.cmt = ifelse.(d.evid .== 1, 1, 3))  # Depot / Central in v2
    read_pumas(d; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
               observations = [:dv], covariates = [:IndividualId])
end

function starts(m, t, s)
    m == "M1" && return M1mod.INITS_FELO
    m == "M2" && return M2_INITS
    m == "M3" && return [l => merge(p, (tvt50 = t, tvgamma = s)) for (l, p) in M3mod.INITS_V2]
    [l => merge(p, (t50R = REF[1], gR = REF[2], t50T = t, gT = s)) for (l, p) in M3mod.INITS_V2]
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

want_m = split(get(ENV, "MC_MODELS", "M1,M2,M3,M4"), ",")
want_a = haskey(ENV, "MC_ARMS") ? Set(split(ENV["MC_ARMS"], ",")) : nothing
for m in want_m, (t, s) in ARMS_ALL
    m == "M4" && (t, s) == REF && continue
    want_a === nothing || akey(t, s) in want_a || continue
    pop = population(m, t, s); model = MODEL[m]
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
                          LL = missing, m2LL = missing, n_par = NPAR[m], n_obs = nobs, AIC = missing, BIC = missing,
                          sigma = missing, secs = secs))
            continue
        end
        ll = loglikelihood(f); k = NPAR[m]
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
