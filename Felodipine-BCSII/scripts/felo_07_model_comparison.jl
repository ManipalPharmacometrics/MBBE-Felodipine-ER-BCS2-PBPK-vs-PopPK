# =============================================================================
# felo_07_model_comparison.jl — gather every felodipine PopPK model fit into one
# comparison folder (2026-10-23). Reads the poppk-fits outputs; no Pumas.
#
# Models (equations in poppk/model_comparison/MODELS.md):
#   M1 1-cmt + Weibull (v1 baseline)        M2 2-cmt + Weibull (10-16, rejected)
#   M3 Weibull -> ka -> 2-cmt, per arm (v2)  M4 M3 paired ref + test, shared disposition
#   M1d/M2d/M3d  M1-M3 in the @delay form, lag 0.3113 h; M1d0.. lag 0 (05_delay_model_fits.jl)
#
# Inputs (FITS = poppk-fits):
#   FITS/output/arms/       M1, all 27 arms (best of 2 starts; LL of the other kept)
#   FITS/output_v2/arms/    M3, all 27 arms (same layout)
#   FITS/output_modelcomp/  M1-M4 on 5 diagnostic arms, every start (04_model_comparison_fits.jl)
#   FITS/output/pbpk_pk_parameters_individual.csv   PK-Sim truth per subject
# Outputs (poppk/model_comparison/):
#   all27_M1_vs_M3.csv        per arm x model: LL, AIC, BIC, sigma, endpoint agreement,
#                             Cmax/AUC GMR vs reference (model) and PBPK GMR
#   diagnostic_fits.csv       every model x arm x start: fit statistics + endpoint agreement
#   diagnostic_params.csv     every model x arm x start: all parameter estimates (long)
#   diagnostic_best.csv       best start per model x arm, wide, for reading
#   delay_vs_ode_fits.csv     M1-M3 vs their @delay versions (best starts): LL, endpoints, GMR
#   delay_vs_ode_params.csv   the same pairs, every parameter side by side with % difference
#   fig_all27_M1_vs_M3.png    Cmax agreement and GMR error by arm
# Usage (from programs):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_07_model_comparison.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie, Printf

const FITS = normpath(joinpath(@__DIR__, "../../poppk-fits"))
const OUT  = joinpath(@__DIR__, "..", "poppk", "model_comparison"); mkpath(OUT)
const REF  = (T50 = 4.78, shape = 0.764)
const NPAR = Dict("M1" => 8, "M2" => 10, "M3" => 11, "M4" => 13)   # θ + Ω diag + σ (lag fixed)
const NOBS = 1800                                                   # 100 subjects x 18 samples per arm

gm(x) = exp(mean(log.(x)))
q(x, p) = exp(quantile(log.(x), p))
pb = CSV.read(joinpath(FITS, "output", "pbpk_pk_parameters_individual.csv"), DataFrame)
pbarm(t, s) = filter(r -> r.T50 == t && r.shape == s, pb)

"Agreement of model endpoints `ep` (IndividualId, Cmax, AUC_tEnd) with PK-Sim for arm (t, s)."
function agreement(ep, t, s)
    j = innerjoin(ep[:, [:IndividualId, :Cmax, :AUC_tEnd]],
                  rename(pbarm(t, s)[:, [:IndividualId, :Cmax, :AUC_tEnd]], :Cmax => :Cp, :AUC_tEnd => :Ap),
                  on = :IndividualId)
    ok = isfinite.(j.Cmax) .& isfinite.(j.AUC_tEnd)       # M4 fits give Inf IPRED for some subjects
    rc = j.Cmax[ok] ./ j.Cp[ok]; ra = j.AUC_tEnd[ok] ./ j.Ap[ok]
    (Cmax_ratio_GM = gm(rc), Cmax_ratio_p05 = q(rc, 0.05), Cmax_ratio_p95 = q(rc, 0.95),
     AUC_ratio_GM = gm(ra), AUC_ratio_p05 = q(ra, 0.05), AUC_ratio_p95 = q(ra, 0.95), n_subj = nrow(j),
     n_nonfinite = count(!, ok))
end

"Paired T/R geometric mean ratio for metric `m` (per subject, as the BE analysis)."
function gmr(test, ref, m)
    j = innerjoin(test[:, [:IndividualId, m]], rename(ref[:, [:IndividualId, m]], m => :r), on = :IndividualId)
    x = j[!, m] ./ j.r
    gm(x[isfinite.(x)])
end

# ---- 1. all 27 arms, M1 vs M3 ------------------------------------------------
function load_perarm(dir, model)
    rows = DataFrame()
    eps = Dict{Tuple{Float64,Float64},DataFrame}()
    for d in sort(readdir(dir))
        isfile(joinpath(dir, d, "nca.csv")) || continue
        p = CSV.read(joinpath(dir, d, "params.csv"), DataFrame)[1, :]
        ep = CSV.read(joinpath(dir, d, "nca.csv"), DataFrame)
        eps[(p.T50, p.shape)] = ep
        k = NPAR[model]; ll = p.LL
        push!(rows, merge((model = model, T50 = p.T50, shape = p.shape, best_start = p.basin, LL = ll,
                           LL_other_start = p.LL_other, AIC = -2ll + 2k, BIC = -2ll + k * log(NOBS),
                           sigma = p.sigma_p), agreement(ep, p.T50, p.shape)); cols = :union)
    end
    ref = get(eps, (REF.T50, REF.shape), nothing)
    pref = pbarm(REF.T50, REF.shape)
    rows.GMR_Cmax_model = [ref === nothing ? missing : gmr(eps[(r.T50, r.shape)], ref, :Cmax) for r in eachrow(rows)]
    rows.GMR_AUC_model  = [ref === nothing ? missing : gmr(eps[(r.T50, r.shape)], ref, :AUC_tEnd) for r in eachrow(rows)]
    rows.GMR_Cmax_PBPK  = [gmr(pbarm(r.T50, r.shape), pref, :Cmax) for r in eachrow(rows)]
    rows.GMR_AUC_PBPK   = [gmr(pbarm(r.T50, r.shape), pref, :AUC_tEnd) for r in eachrow(rows)]
    rows.GMR_Cmax_err_pct = 100 .* (rows.GMR_Cmax_model ./ rows.GMR_Cmax_PBPK .- 1)
    rows
end

a27 = vcat(load_perarm(joinpath(FITS, "output", "arms"), "M1"),
           isdir(joinpath(FITS, "output_v2", "arms")) ? load_perarm(joinpath(FITS, "output_v2", "arms"), "M3") : DataFrame();
           cols = :union)
sort!(a27, [:shape, :T50, :model])
CSV.write(joinpath(OUT, "all27_M1_vs_M3.csv"), a27)

# ---- 2. diagnostic arms, every model x start ---------------------------------
const MC = joinpath(FITS, "output_modelcomp")
fits = DataFrame(); params = DataFrame(); refeps = Dict{Tuple{String,String},DataFrame}()
if isdir(MC)
    for m in sort(readdir(MC)), a in sort(readdir(joinpath(MC, m))), st in sort(readdir(joinpath(MC, m, a)))
        dir = joinpath(MC, m, a, st)
        isfile(joinpath(dir, "fitstats.csv")) || continue
        fs = CSV.read(joinpath(dir, "fitstats.csv"), DataFrame)[1, :]
        row = Dict{Symbol,Any}(pairs(NamedTuple(fs)))
        if isfile(joinpath(dir, "endpoints.csv"))
            ep = CSV.read(joinpath(dir, "endpoints.csv"), DataFrame; types = Dict(:period => String))   # "T" would parse as Bool
            for per in unique(ep.period)
                e = filter(r -> r.period == per, ep)
                ag = agreement(e, e.T50[1], e.shape[1])
                for (k, v) in pairs(ag); row[Symbol(k, "_", per)] = v; end
            end
            m == "M4" && (row[:GMR_Cmax_model] = gmr(filter(r -> r.period == "T", ep), filter(r -> r.period == "R", ep), :Cmax))
            m != "M4" && a == "$(REF.T50)_$(REF.shape)" && (refeps[(m, st)] = ep)
        end
        if isfile(joinpath(dir, "coef.csv"))
            c = CSV.read(joinpath(dir, "coef.csv"), DataFrame)
            insertcols!(c, 1, :model => m, :arm => a, :start => st)
            append!(params, c; cols = :union)
        end
        push!(fits, NamedTuple(row); cols = :union)
    end
end
if nrow(fits) > 0
    # per-arm models: T/R GMR against the same model's best reference fit
    hasproperty(fits, :GMR_Cmax_model) || (fits.GMR_Cmax_model = Vector{Union{Missing,Float64}}(missing, nrow(fits)))
    for r in eachrow(fits)
        r.model == "M4" && continue
        isfile(joinpath(MC, r.model, r.arm, r.start, "endpoints.csv")) || continue   # skipped/failed fits
        cands = [(st, ep) for ((m, st), ep) in refeps if m == r.model]
        isempty(cands) && continue
        bestst = fits[(fits.model .== r.model) .& (fits.arm .== "$(REF.T50)_$(REF.shape)"), :]
        bs = bestst.start[argmax(coalesce.(bestst.LL, -Inf))]
        ep = CSV.read(joinpath(MC, r.model, r.arm, r.start, "endpoints.csv"), DataFrame; types = Dict(:period => String))
        r.GMR_Cmax_model = gmr(ep, refeps[(r.model, bs)], :Cmax)
    end
    t50 = [parse(Float64, split(a, "_")[1]) for a in fits.arm]; sh = [parse(Float64, split(a, "_")[2]) for a in fits.arm]
    fits.GMR_Cmax_PBPK = [gmr(pbarm(t, s), pbarm(REF.T50, REF.shape), :Cmax) for (t, s) in zip(t50, sh)]
    fits.GMR_Cmax_err_pct = 100 .* (coalesce.(fits.GMR_Cmax_model, NaN) ./ fits.GMR_Cmax_PBPK .- 1)
    # ΔAIC within arm (M4 is fitted on two arms of data, so not comparable: excluded)
    fits.dAIC_vs_best_in_arm = [r.model == "M4" ? missing :
        r.AIC - minimum(skipmissing(fits[(fits.model .!= "M4") .& (fits.arm .== r.arm), :AIC]); init = Inf)
        for r in eachrow(fits)]
    sort!(fits, [:arm, :model, :start])
    CSV.write(joinpath(OUT, "diagnostic_fits.csv"), fits)
    CSV.write(joinpath(OUT, "diagnostic_params.csv"), params)
    # best start per model x arm, parameters wide
    best = combine(groupby(fits, [:model, :arm]), sdf -> sdf[argmax(coalesce.(sdf.LL, -Inf)), :])
    wide = unstack(params, [:model, :arm, :start], :param, :value)
    bw = leftjoin(best, wide, on = [:model, :arm, :start])
    sort!(bw, [:arm, :model])
    CSV.write(joinpath(OUT, "diagnostic_best.csv"), bw)

    # ODE hazard (M1-M3) vs @delay form (M1d-M3d): best start of each, per arm
    dvo = DataFrame(); dvo_par = DataFrame()
    for base in ("M1", "M2", "M3"), v in ("d", "d0"), a in unique(best.arm)
        o = best[(best.model .== base) .& (best.arm .== a), :]; d = best[(best.model .== base * v) .& (best.arm .== a), :]
        (nrow(o) == 1 && nrow(d) == 1) || continue
        o = o[1, :]; d = d[1, :]
        push!(dvo, (model = base, delay_variant = base * v, lag_delay = v == "d" ? 0.3113 : 0.0, arm = a, start_ode = o.start, start_delay = d.start,
                    status_ode = o.status, status_delay = d.status, LL_ode = o.LL, LL_delay = d.LL,
                    dLL = d.LL - o.LL, sigma_ode = o.sigma, sigma_delay = d.sigma,
                    Cmax_ratio_ode = o.Cmax_ratio_GM_T, Cmax_ratio_delay = d.Cmax_ratio_GM_T,
                    AUC_ratio_ode = o.AUC_ratio_GM_T, AUC_ratio_delay = d.AUC_ratio_GM_T,
                    GMR_Cmax_ode = o.GMR_Cmax_model, GMR_Cmax_delay = d.GMR_Cmax_model, GMR_Cmax_PBPK = o.GMR_Cmax_PBPK,
                    secs_ode = o.secs, secs_delay = d.secs))
        po = params[(params.model .== base) .& (params.arm .== a) .& (params.start .== o.start), :]
        pd = params[(params.model .== base * v) .& (params.arm .== a) .& (params.start .== d.start), :]
        j = outerjoin(rename(po[:, [:param, :value]], :value => :ode), rename(pd[:, [:param, :value]], :value => :delay), on = :param)
        j.pct_diff = 100 .* (j.delay ./ j.ode .- 1)
        insertcols!(j, 1, :model => base, :delay_variant => base * v, :arm => a)
        append!(dvo_par, j; cols = :union)
    end
    if nrow(dvo) > 0
        CSV.write(joinpath(OUT, "delay_vs_ode_fits.csv"), sort!(dvo, [:model, :delay_variant, :arm]))
        CSV.write(joinpath(OUT, "delay_vs_ode_params.csv"), sort!(dvo_par, [:model, :delay_variant, :arm, :param]))
    end
end

# ---- 3. figure: all 27 arms ----------------------------------------------------
if "M3" in a27.model
    COL = Dict("M1" => colorant"#eb6834", "M3" => colorant"#2a78d6")
    lab = Dict("M1" => "M1  1-cmt + Weibull (v1)", "M3" => "M3  Weibull → ka → 2-cmt (v2)")
    arms = sort(unique(a27[:, [:T50, :shape]]), [:shape, :T50])
    xl = ["$(r.T50)/$(r.shape)" for r in eachrow(arms)]
    idx(r) = findfirst(i -> arms.T50[i] == r.T50 && arms.shape[i] == r.shape, 1:nrow(arms))
    fig = Figure(size = (1500, 820), fontsize = 14)
    ax1 = Axis(fig[1, 1], ylabel = "Cmax PopPK / PBPK (GM, p05–p95)", xticks = (1:nrow(arms), xl),
               xticklabelrotation = π / 3, title = "Per-subject Cmax agreement (all 27 arms)")
    ax2 = Axis(fig[2, 1], ylabel = "Cmax GMR error vs PBPK (%)", xticks = (1:nrow(arms), xl),
               xticklabelrotation = π / 3, xlabel = "arm (T50 h / shape)", title = "Test/reference Cmax GMR: PopPK vs PBPK")
    hlines!(ax1, [1.0], color = :gray50); hlines!(ax2, [0.0], color = :gray50)
    for (k, m) in enumerate(["M1", "M3"])
        d = a27[a27.model .== m, :]; x = idx.(eachrow(d)) .+ (k == 1 ? -0.15 : 0.15)
        rangebars!(ax1, x, d.Cmax_ratio_p05, d.Cmax_ratio_p95, color = COL[m])
        scatter!(ax1, x, d.Cmax_ratio_GM, color = COL[m], label = lab[m])
        scatter!(ax2, x, d.GMR_Cmax_err_pct, color = COL[m], label = lab[m])
    end
    axislegend(ax1, position = :lt);
    save(joinpath(OUT, "fig_all27_M1_vs_M3.png"), fig, px_per_unit = 1.5)
end

# ---- console summary -----------------------------------------------------------
for m in unique(a27.model)
    d = a27[a27.model .== m, :]
    @printf("%s  arms %2d  ΣLL %9.1f  Cmax ratio %.3f–%.3f  |GMR err| max %.1f%%  mean %.1f%%\n", m, nrow(d),
            sum(d.LL), minimum(d.Cmax_ratio_GM), maximum(d.Cmax_ratio_GM),
            maximum(abs.(skipmissing(d.GMR_Cmax_err_pct))), mean(abs.(skipmissing(d.GMR_Cmax_err_pct))))
end
nrow(fits) > 0 && show(select(fits, :arm, :model, :start, :status, :LL, :AIC, :sigma,
                              r"Cmax_ratio_GM", :GMR_Cmax_err_pct), allrows = true, allcols = true)
println("\nWrote -> $OUT")
