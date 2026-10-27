# =============================================================================
# 02_fit_poppk.jl — fit the felodipine PopPK model to each formulation arm
# independently (FOCE), then NCA on the dense individual predictions (IPRED).
# Mirrors the metformin per-arm 02 + 03.
#
# Resumable: one folder per arm under output/arms/, skipped if done.
# Optional env: POPPK_ARMS = "4.78_0.764,2.39_0.764" to fit only those arms;
#               POPPK_SHARD = "k/m" to fit every m-th arm starting at k (parallel).
# Per arm writes: params.csv, se.csv, shrinkage.csv, nca.csv (per subject
#   Cmax nM, AUC_tEnd nM·h over 0-96 h, Tmax), ipred_thin.csv (0.5-h grid).
# Run: julia -J/opt/julia/lib/julia/pumas.so --project=<repo>/julia-env scripts/02_fit_poppk.jl
# =============================================================================

using Pumas, CSV, DataFrames, Statistics

const HERE = @__DIR__
const ROOT = normpath(joinpath(HERE, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output"))
const ARMS = joinpath(OUT, "arms"); mkpath(ARMS)
include(joinpath(ROOT, "src", "poppk_model_felo_1cmt.jl"))   # baseline (2-cmt not identifiable, notes §9.26)

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

trap(t, c) = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in 2:length(t))

for a in eachrow(arms)
    dir = joinpath(ARMS, key(a.T50, a.shape))
    isfile(joinpath(dir, "nca.csv")) && (println("skip ", key(a.T50, a.shape)); continue)
    mkpath(dir)
    arm = filter(r -> r.T50 == a.T50 && r.shape == a.shape, data)
    pop = read_pumas(arm; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                     observations = [:dv], covariates = [:T50, :shape, :IndividualId])
    @info "Fitting T50 = $(a.T50) h, shape = $(a.shape) (n = $(length(pop)))"
    t0 = time()
    # two starting points (release- vs disposition-limited basin); keep higher LL
    fits = []
    for (lab, p) in INITS_FELO
        f = try
            fit(poppk_model_felo, pop, merge(init_params(poppk_model_felo), p), FOCE();
                constantcoef = (:tvlag,), optim_options = (show_trace = false,))
        catch err
            @warn "fit failed" lab err; nothing
        end
        f === nothing || push!(fits, (lab = lab, f = f, ll = loglikelihood(f)))
    end
    best = fits[argmax([x.ll for x in fits])]
    fpm = best.f
    el = round(time() - t0, digits = 1)
    c = coef(fpm)
    ofv = try; deviance(fpm); catch; missing; end
    CSV.write(joinpath(dir, "params.csv"),
        DataFrame(T50 = a.T50, shape = a.shape, tvalpha = c.tvalpha, tvgamma = c.tvgamma,
                  T50_implied = c.tvalpha * log(2)^(1 / c.tvgamma),
                  tvcl = c.tvcl, tvvc = c.tvvc, tvlag = c.tvlag,
                  omega_alpha = c.Ω[1, 1], omega_cl = c.Ω[2, 2], omega_vc = c.Ω[3, 3],
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
        DataFrame(T50 = a.T50, shape = a.shape, eta_alpha = ηs[1], eta_cl = ηs[2], eta_vc = ηs[3]))
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
    println("done $(key(a.T50, a.shape)) [$(best.lab)] in $(el) s  T50_implied=$(round(c.tvalpha*log(2)^(1/c.tvgamma), digits=2))  gamma=$(round(c.tvgamma, digits=3))  CL=$(round(c.tvcl, digits=1))  sigma=$(round(c.σ_p, digits=3))")
end
