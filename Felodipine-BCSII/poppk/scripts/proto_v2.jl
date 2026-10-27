# proto_v2.jl — prototype check of the v2 model on diagnostic arms (per-arm).
# Reports fitted T50/shape vs PK-Sim truth and IPRED/PBPK Cmax, AUC GM.
# Env: PROTO_ARMS = "4.78_0.764,4.78_0.5,4.78_1.5" (default)
using Pumas, CSV, DataFrames, Statistics
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = joinpath(ROOT, "output"); const PDIR = joinpath(OUT, "proto_v2"); mkpath(PDIR)
include(joinpath(ROOT, "src", "poppk_model_felo_v2.jl"))
data = CSV.read(joinpath(OUT, "be_dataset.csv"), DataFrame)
pb   = CSV.read(joinpath(OUT, "pbpk_pk_parameters_individual.csv"), DataFrame)
want = split(get(ENV, "PROTO_ARMS", "4.78_0.764,4.78_0.5,4.78_1.5"), ",")
const DENSE = collect(0.0:0.05:96.0)
trap(t, c) = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in 2:length(t))
gm(x) = exp(mean(log.(x)))
res = DataFrame()
for w in want
    T, S = parse.(Float64, split(w, "_"))
    arm = filter(r -> r.T50 == T && r.shape == S, data)
    arm.cmt = ifelse.(arm.evid .== 1, 1, 3)    # Depot / Central index in v2 (Depot,Gut,Central,Periph)
    pop = read_pumas(arm; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                     observations = [:dv], covariates = [:IndividualId])
    for (lab, p) in INITS_V2
        t0 = time()
        f = try
            fit(poppk_felo_wka2, pop, merge(init_params(poppk_felo_wka2), p, (tvt50 = T, tvgamma = S)), FOCE();
                constantcoef = (:tvlag,), optim_options = (show_trace = false,))
        catch err
            @warn "fit failed" w lab err; nothing
        end
        f === nothing && continue
        c = coef(f); ll = loglikelihood(f)
        pr = predict(f; obstimes = DENSE)
        e = DataFrame(IndividualId = [Int(s.covariates(0.0).IndividualId) for s in pop],
                      Cmax = [maximum(x.ipred.dv) for x in pr],
                      AUC = [trap(x.time, x.ipred.dv) for x in pr])
        j = innerjoin(e, filter(r -> r.T50 == T && r.shape == S, pb), on = :IndividualId, makeunique = true)
        r = (arm = w, start = lab, LL = ll, secs = round(time() - t0), t50 = c.tvt50, gamma = c.tvgamma,
             ka = c.tvka, cl = c.tvcl, vc = c.tvvc, q = c.tvq, vp = c.tvvp, sigma = c.σ_p,
             om = round.(diag(c.Ω), digits = 3),
             Cmax_ratio = gm(j.Cmax ./ j.Cmax_1), AUC_ratio = gm(j.AUC ./ j.AUC_tEnd))
        push!(res, r; cols = :union); println(r)
        CSV.write(joinpath(PDIR, "proto_wka2.csv"), res)
    end
end
