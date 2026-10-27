# proto_pair.jl — prototype of the paired (ref + test, same subjects) v2 fit.
# Test period dosed at TP2 into Depot2. Reports per-period IPRED/PBPK Cmax, AUC
# and the T/R Cmax GMR from PopPK vs PBPK.
# Env: PROTO_ARMS = "4.78_0.5,4.78_1.5" (test arms; reference is 4.78_0.764)
using Pumas, CSV, DataFrames, Statistics
const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = joinpath(ROOT, "output"); const PDIR = joinpath(OUT, "proto_v2"); mkpath(PDIR)
include(joinpath(ROOT, "src", "poppk_model_felo_v2.jl"))
data = CSV.read(joinpath(OUT, "be_dataset.csv"), DataFrame)
pb   = CSV.read(joinpath(OUT, "pbpk_pk_parameters_individual.csv"), DataFrame)
const REF = (4.78, 0.764)
want = split(get(ENV, "PROTO_ARMS", "4.78_0.5,4.78_1.5"), ",")
const DENSE = vcat(collect(0.0:0.05:96.0), collect(TP2:0.05:TP2 + 96.0))
trap(t, c) = sum((t[i] - t[i-1]) * (c[i] + c[i-1]) / 2 for i in 2:length(t))
gm(x) = exp(mean(log.(x)))

function pair_data(T, S)
    r = filter(x -> x.T50 == REF[1] && x.shape == REF[2], data)
    t = filter(x -> x.T50 == T && x.shape == S, data)
    r.cmt = ifelse.(r.evid .== 1, 1, 4)
    t.cmt = ifelse.(t.evid .== 1, 2, 4); t.time .+= TP2
    d = vcat(r, t); d.ID = d.IndividualId; sort!(d, [:ID, :time]); d
end

res = DataFrame()
for w in want
    T, S = parse.(Float64, split(w, "_"))
    pop = read_pumas(pair_data(T, S); id = :ID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                     observations = [:dv])
    for (lab, p) in INITS_V2
        t0 = time()
        f = try
            fit(poppk_felo_pair, pop,
                merge(init_params(poppk_felo_pair), p, (t50R = REF[1], gR = REF[2], t50T = T, gT = S)), FOCE();
                constantcoef = (:tvlag,), optim_options = (show_trace = false,))
        catch err
            @warn "fit failed" w lab err; nothing
        end
        f === nothing && continue
        c = coef(f); ll = loglikelihood(f)
        pr = predict(f; obstimes = DENSE)
        e = DataFrame(IndividualId = Int[], CmaxR = Float64[], CmaxT = Float64[], AUCR = Float64[], AUCT = Float64[])
        for (s, x) in zip(pop, pr)
            tt = x.time; ip = x.ipred.dv; a = tt .<= 96.0; b = tt .>= TP2
            push!(e, (parse(Int, string(s.id)), maximum(ip[a]), maximum(ip[b]), trap(tt[a], ip[a]), trap(tt[b] .- TP2, ip[b])))
        end
        pr_ = filter(r -> r.T50 == REF[1] && r.shape == REF[2], pb)
        pt_ = filter(r -> r.T50 == T && r.shape == S, pb)
        j = innerjoin(innerjoin(e, rename(pr_[:, [:IndividualId, :Cmax, :AUC_tEnd]], :Cmax => :pCR, :AUC_tEnd => :pAR), on = :IndividualId),
                      rename(pt_[:, [:IndividualId, :Cmax, :AUC_tEnd]], :Cmax => :pCT, :AUC_tEnd => :pAT), on = :IndividualId)
        r = (arm = w, start = lab, LL = ll, secs = round(time() - t0), t50R = c.t50R, gR = c.gR, t50T = c.t50T, gT = c.gT,
             ka = c.tvka, cl = c.tvcl, vc = c.tvvc, q = c.tvq, vp = c.tvvp, sigma = c.σ_p,
             CmaxR_ratio = gm(j.CmaxR ./ j.pCR), CmaxT_ratio = gm(j.CmaxT ./ j.pCT),
             AUCR_ratio = gm(j.AUCR ./ j.pAR), AUCT_ratio = gm(j.AUCT ./ j.pAT),
             GMR_Cmax_popPK = gm(j.CmaxT ./ j.CmaxR), GMR_Cmax_PBPK = gm(j.pCT ./ j.pCR))
        push!(res, r; cols = :union); println(r)
        CSV.write(joinpath(PDIR, "proto_pair.csv"), res)
    end
end
