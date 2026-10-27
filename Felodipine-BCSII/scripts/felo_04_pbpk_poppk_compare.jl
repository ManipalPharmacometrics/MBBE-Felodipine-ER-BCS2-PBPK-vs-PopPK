# =============================================================================
# felo_04_pbpk_poppk_compare.jl
# Felodipine ER (MBBE drug two): PBPK vs PopPK decision equivalence along the
# two CQA axes, plus Type I error. Felodipine version of the metformin
# 04_be_safespace_overlay.jl + 05_type1_overlay.jl, using the SAME BEUtils.
#
# Axes (1-D slices of the 27-formulation design):
#   T50   axis: shape = 0.764, T50 2.39-9.56 h (9 anchors; log-linear interp.)
#   shape axis: T50 = 4.78 h, shape 0.5-1.5   (7 anchors; linear in shape)
# Endpoints: Cmax, AUC_0-t (AUC_tEnd, 0-96 h), joint (both). No AUC_inf
# (publication convention). WSV: felodipine fasted (Cmax 35.7%, AUC_0-t 25.3%).
# n = 24 and 54 (54 = powered 2x2 n at GMR 0.95).
#
# Type I error (Morris): invert CQA -> GMR to find the simulated BE boundary
# (first exit from [0.80, 1.25] walking out from the reference), then the
# Monte-Carlo acceptance rate there (nominal 5%). Cross-evaluation: each
# method's acceptance at the OTHER method's boundary. Plus the TOST
# calibration control (reference x 0.80) and, for the HVD story, the
# acceptance of FDA RSABE / EMA ABEL (n = 24, Cmax) at the ABE boundary.
#
# Inputs : <dir>/pbpk_pk_parameters_individual.csv, <dir>/poppk_pk_parameters_individual.csv
# Outputs: <dir>/compare/*.csv and fig_*.png
# Usage (from programs):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_04_pbpk_poppk_compare.jl ../Felodipine-BCSII/poppk/output_rich
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "..", "programs", "julia", "be_utils.jl"))
using .BEUtils

const DIR  = ARGS[1]
const OUT  = joinpath(DIR, "compare"); mkpath(OUT)
const CV   = wsv_cv(:felodipine, :fasted)
const CVWR = wsv_cv(:felodipine, :cvwr_rsabe, "Cmax")
const NS   = [24, 54]
const N_T1 = parse(Int, get(ENV, "N_T1", "20000"))   # Type I: SE ~0.15 pp at 5%
const N_PW = parse(Int, get(ENV, "N_PW", "1000"))    # power / OC curves
const REF  = (T50 = 4.78, shape = 0.764)
const SEED = 1234                                     # VBE resampling seed (be_utils default, stated explicitly)

pb = CSV.read(joinpath(DIR, "pbpk_pk_parameters_individual.csv"), DataFrame)
pp = CSV.read(joinpath(DIR, "poppk_pk_parameters_individual.csv"), DataFrame)
SOURCES = ["PBPK" => pb, "PopPK" => pp]
COL = Dict("PBPK" => colorant"#2a78d6", "PopPK" => colorant"#eb6834")
INK, INK2 = colorant"#1f1f1d", colorant"#5f5e5a"

# ---- axes ----------------------------------------------------------------------
AXES = [(name = "T50", col = :T50, fix = :shape, fixv = REF.shape, ref = REF.T50,
         logx = true, label = "Dissolution T50 (h), shape 0.764"),
        (name = "shape", col = :shape, fix = :T50, fixv = REF.T50, ref = REF.shape,
         logx = false, label = "Weibull shape, T50 4.78 h")]
slice(pk, ax) = filter(r -> r[ax.fix] == ax.fixv, pk)
anchors(pk, ax) = sort(unique(slice(pk, ax)[!, ax.col]))
ids_of(pk) = sort(unique(pk.IndividualId))

function metric_matrix(pk, ax, metric)
    s = slice(pk, ax); A = anchors(pk, ax); I = ids_of(pk)
    M = Matrix{Float64}(undef, length(I), length(A))
    for (j, a) in enumerate(A)
        d = sort(filter(r -> r[ax.col] == a, s), :IndividualId)
        @assert d.IndividualId == I
        M[:, j] = d[!, metric]
    end
    M, A
end
function interp(vals, A, x, logx)
    x <= A[1] && return vals[1]; x >= A[end] && return vals[end]
    k = searchsortedlast(A, x); x == A[k] && return vals[k]
    w = logx ? (log(x) - log(A[k])) / (log(A[k+1]) - log(A[k])) : (x - A[k]) / (A[k+1] - A[k])
    exp((1 - w) * log(vals[k]) + w * log(vals[k+1]))
end
function test_ref_at(pk, ax, metric, x)
    M, A = metric_matrix(pk, ax, metric)
    rj = findfirst(==(ax.ref), A)
    ([interp(M[i, :], A, x, ax.logx) for i in axes(M, 1)], M[:, rj])
end
gmr_at(pk, ax, metric, x) = be_paired(test_ref_at(pk, ax, metric, x)...).GMR

function accept(pk, ax, ep, x, n; n_trial = N_T1)
    if ep == "Joint"
        tC, rC = test_ref_at(pk, ax, "Cmax", x); tA, rA = test_ref_at(pk, ax, "AUC_tEnd", x)
        vbe_replicate_joint(tC, rC, tA, rA; seed = SEED, n_subj = n, n_trial, cv_C = CV["Cmax"],
                            cv_A = CV["AUC_tEnd"]).pass_rate
    else
        t, r = test_ref_at(pk, ax, ep, x)
        vbe_replicate(t, r; seed = SEED, n_subj = n, n_trial, wsv_cv = CV[ep]).pass_rate
    end
end
function accept_scaled(pk, ax, x, method; n = 24, n_trial = N_T1)
    t, r = test_ref_at(pk, ax, "Cmax", x)
    vbe_replicate_scaled(t, r; seed = SEED, n_subj = n, n_trial, cv_wr = CVWR, method).pass_rate
end
function wilson(p, n; z = 1.96)
    d = 1 + z^2 / n; c = (p + z^2 / (2n)) / d
    h = z * sqrt(p * (1 - p) / n + z^2 / (4n^2)) / d
    (max(0.0, c - h), min(1.0, c + h))
end
grid_of(A, logx) = logx ? exp.(range(log(A[1]), log(A[end]), length = 161)) :
                          collect(range(A[1], A[end], length = 161))

"First CQA value, walking out from the reference, where the GMR crosses 0.80 or 1.25."
function first_exit(G, g, ref, side)
    ri = argmin(abs.(G .- ref)); step = side === :low ? -1 : 1; i = ri
    while side === :low ? i > 1 : i < length(G)
        j = i + step
        for (nm, tg) in (("0.80", BE_LOWER), ("1.25", BE_UPPER))
            a, b = g[i], g[j]
            if (a - tg) * (b - tg) <= 0 && a != b
                w = (tg - a) / (b - a)
                return (x = G[i] + w * (G[j] - G[i]), limit = nm)
            end
        end
        i = j
    end
    (x = missing, limit = "not reached")
end

# ---- 1. GMR + 90% CI at the anchors -----------------------------------------------
gmr = DataFrame()
for (src, pk) in SOURCES, ax in AXES, m in ("Cmax", "AUC_tEnd")
    M, A = metric_matrix(pk, ax, m); rj = findfirst(==(ax.ref), A)
    for (j, a) in enumerate(A)
        r = be_paired(M[:, j], M[:, rj])
        push!(gmr, (source = src, axis = ax.name, x = a, metric = m, GMR = r.GMR,
                    CI_lower = r.CI_lower, CI_upper = r.CI_upper))
    end
end
CSV.write(joinpath(OUT, "gmr_axes.csv"), gmr)

# ---- 2. simulated BE boundaries (CQA inversion) -------------------------------------
bnd = DataFrame()
for ax in AXES, (src, pk) in SOURCES
    A = anchors(pk, ax); G = grid_of(A, ax.logx)
    gC = [gmr_at(pk, ax, "Cmax", x) for x in G]; gA = [gmr_at(pk, ax, "AUC_tEnd", x) for x in G]
    for side in (:low, :high)
        eC = first_exit(G, gC, ax.ref, side); eA = first_exit(G, gA, ax.ref, side)
        for (ep, e, who) in (("Cmax", eC, "Cmax"), ("AUC_tEnd", eA, "AUC_tEnd"))
            push!(bnd, (axis = ax.name, source = src, endpoint = ep, side = String(side),
                        bnd = e.x, limit = e.limit, binding = who); cols = :union)
        end
        c = filter(t -> !ismissing(t[1].x), [(eC, "Cmax"), (eA, "AUC_tEnd")])
        (e, who) = isempty(c) ? ((x = missing, limit = "not reached"), "none") :
                   argmin(t -> abs(t[1].x - ax.ref), c)
        push!(bnd, (axis = ax.name, source = src, endpoint = "Joint", side = String(side),
                    bnd = e.x, limit = e.limit, binding = who); cols = :union)
    end
end
CSV.write(joinpath(OUT, "boundaries.csv"), bnd)
println("\nSimulated BE boundaries (GMR leaves [0.80, 1.25]):")
show(unstack(select(bnd, :axis, :endpoint, :side, :source, :bnd), :source, :bnd); allrows = true); println()

# ---- 3. Type I error: own boundary + cross-evaluation ------------------------------
t1 = DataFrame()
for ax in AXES, ep in ("Cmax", "AUC_tEnd", "Joint"), side in ("low", "high"), n in NS
    rows = filter(r -> r.axis == ax.name && r.endpoint == ep && r.side == side, bnd)
    for (src, pk) in SOURCES
        me = only(filter(r -> r.source == src, rows)); ismissing(me.bnd) && continue
        oth = only(filter(r -> r.source != src, rows))
        own = accept(pk, ax, ep, me.bnd, n); lo, hi = wilson(own, N_T1)
        cross = ismissing(oth.bnd) ? missing : accept(pk, ax, ep, oth.bnd, n)
        gx = ismissing(oth.bnd) ? missing :
             gmr_at(pk, ax, me.binding == "none" ? "Cmax" : me.binding, oth.bnd)
        push!(t1, (axis = ax.name, endpoint = ep, side = side, n = n, judged_by = src,
                   own_bnd = me.bnd, limit = me.limit, binding = me.binding,
                   T1_own = own, T1_lo = lo, T1_hi = hi, other_bnd = oth.bnd,
                   T1_at_other_bnd = cross, GMR_at_other_bnd = gx); cols = :union)
    end
end
CSV.write(joinpath(OUT, "type1.csv"), t1)
println("\nType I error ($N_T1 trials), PBPK truth at PopPK's boundary = T1_at_other_bnd on PBPK rows:")
show(select(t1, :axis, :endpoint, :side, :n, :judged_by, :own_bnd, :T1_own, :other_bnd,
            :T1_at_other_bnd, :GMR_at_other_bnd); allrows = true); println()

# ---- 4. calibration control + scaled methods at the ABE Cmax boundary ---------------
ctl = DataFrame()
for (src, pk) in SOURCES, ep in ("Cmax", "AUC_tEnd", "Joint"), n in NS
    ax = AXES[1]
    p = if ep == "Joint"
        _, rC = test_ref_at(pk, ax, "Cmax", ax.ref); _, rA = test_ref_at(pk, ax, "AUC_tEnd", ax.ref)
        vbe_replicate_joint(rC .* BE_LOWER, rC, rA .* BE_LOWER, rA; seed = SEED, n_subj = n, n_trial = N_T1,
                            cv_C = CV["Cmax"], cv_A = CV["AUC_tEnd"]).pass_rate
    else
        _, r = test_ref_at(pk, ax, ep, ax.ref)
        vbe_replicate(r .* BE_LOWER, r; seed = SEED, n_subj = n, n_trial = N_T1, wsv_cv = CV[ep]).pass_rate
    end
    push!(ctl, (source = src, endpoint = ep, n = n, T1_scaled_ref = p))
end
CSV.write(joinpath(OUT, "type1_calibration_control.csv"), ctl)

sc = DataFrame()
for ax in AXES, (src, pk) in SOURCES, side in ("low", "high")
    me = only(filter(r -> r.axis == ax.name && r.source == src && r.endpoint == "Cmax" &&
                          r.side == side, bnd))
    ismissing(me.bnd) && continue
    for method in (:fda_rsabe, :ema_abel)
        push!(sc, (axis = ax.name, source = src, side = side, bnd = me.bnd, limit = me.limit,
                   method = String(method), n = 24, accept = accept_scaled(pk, ax, me.bnd, method)))
    end
end
CSV.write(joinpath(OUT, "scaled_accept_at_abe_boundary.csv"), sc)
println("\nRSABE / ABEL acceptance at the ABE Cmax boundary (true GMR 0.80 or 1.25):")
show(sc; allrows = true); println()

# ---- 5. figures ----------------------------------------------------------------------
set_theme!(fontsize = 15, Axis = (xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
                                  topspinevisible = false, rightspinevisible = false))
ticks_of(ax) = ax.logx ? ([2.39, 3.35, 4.78, 7.17, 9.56], ["2.39", "3.35", "4.78", "7.17", "9.56"]) :
                          ([0.5, 0.764, 1.0, 1.2, 1.5], ["0.5", "0.764", "1.0", "1.2", "1.5"])
dodge(src, ax) = (src == "PBPK" ? -1 : 1) * (ax.logx ? 0.012 : 0.008)

# (a) GMR vs CQA, both axes x {Cmax, AUC_0-t}
let fig = Figure(size = (1200, 760))
    Label(fig[0, 1:2], "Test/reference GMR along the two CQA axes — PBPK vs PopPK (90% CI, 100 paired subjects)";
          font = :bold, fontsize = 17, halign = :left)
    for (i, m) in enumerate(("Cmax", "AUC_tEnd")), (j, ax) in enumerate(AXES)
        a = Axis(fig[i, j]; xlabel = ax.label, ylabel = (m == "Cmax" ? "Cmax" : "AUC₀₋ₜ") * " GMR",
                 xscale = ax.logx ? log10 : identity,
                 xticks = ticks_of(ax))
        hspan!(a, BE_LOWER, BE_UPPER; color = (:black, 0.05))
        hlines!(a, [BE_LOWER, BE_UPPER]; color = INK2, linestyle = :dash, linewidth = 1)
        vlines!(a, [ax.ref]; color = INK2, linestyle = :dot, linewidth = 1)
        for (src, _) in SOURCES
            d = sort(filter(r -> r.source == src && r.axis == ax.name && r.metric == m, gmr), :x)
            xs = ax.logx ? d.x .* exp(dodge(src, ax)) : d.x .+ dodge(src, ax)
            errorbars!(a, xs, d.GMR, d.GMR .- d.CI_lower, d.CI_upper .- d.GMR;
                       color = COL[src], whiskerwidth = 6)
            scatterlines!(a, xs, d.GMR; color = COL[src], markersize = 9, linewidth = 2, label = src)
        end
    end
    Legend(fig[3, 1:2], [LineElement(color = COL[s], linewidth = 3) for s in first.(SOURCES)],
           first.(SOURCES); orientation = :horizontal, framevisible = false)
    save(joinpath(OUT, "fig_gmr_axes.png"), fig, px_per_unit = 2)
end

# (b) joint ABE power along both axes, n = 24 / 54
pw = DataFrame()
for ax in AXES, (src, pk) in SOURCES, n in NS
    G = grid_of(anchors(pk, ax), ax.logx)[1:2:end]
    for x in G
        push!(pw, (axis = ax.name, source = src, n = n, x = x,
                   power = accept(pk, ax, "Joint", x, n; n_trial = N_PW)))
    end
end
CSV.write(joinpath(OUT, "power_joint_axes.csv"), pw)
let fig = Figure(size = (1200, 760))
    Label(fig[0, 1:2], "Joint ABE power (Cmax and AUC₀₋ₜ, 2×2) along the CQA axes — PBPK vs PopPK ($N_PW trials per point)";
          font = :bold, fontsize = 17, halign = :left)
    for (i, n) in enumerate(reverse(NS)), (j, ax) in enumerate(AXES)
        a = Axis(fig[i, j]; xlabel = ax.label, ylabel = "P(conclude BE), n = $n",
                 xscale = ax.logx ? log10 : identity,
                 xticks = ticks_of(ax))
        ylims!(a, -0.02, 1.02)
        hlines!(a, [0.05, 0.80]; color = INK2, linestyle = :dot, linewidth = 1)
        vlines!(a, [ax.ref]; color = INK2, linestyle = :dot, linewidth = 1)
        for (src, _) in SOURCES
            d = sort(filter(r -> r.source == src && r.axis == ax.name && r.n == n, pw), :x)
            lines!(a, d.x, d.power; color = COL[src], linewidth = 2.5, label = src)
        end
    end
    Legend(fig[3, 1:2], [LineElement(color = COL[s], linewidth = 3) for s in first.(SOURCES)],
           first.(SOURCES); orientation = :horizontal, framevisible = false)
    save(joinpath(OUT, "fig_power_joint_axes.png"), fig, px_per_unit = 2)
end

# (c) operating characteristic along the T50 axis: acceptance vs the true Cmax GMR
#     (the shape axis never reaches GMR 0.80, so it adds nothing here)
let fig = Figure(size = (1200, 500)), ax = AXES[1]
    Label(fig[0, 1:2], "Operating characteristic — P(conclude BE, joint ABE) vs the true Cmax GMR, T50 axis";
          font = :bold, fontsize = 17, halign = :left)
    for (j, n) in enumerate(NS)
        a = Axis(fig[1, j]; xlabel = "true Cmax GMR (test/reference)",
                 ylabel = j == 1 ? "P(conclude BE)" : "", title = "n = $n", titlealign = :left)
        ylims!(a, -0.02, 1.02)
        vlines!(a, [BE_LOWER, BE_UPPER]; color = INK2, linestyle = :dash, linewidth = 1)
        hlines!(a, [0.05]; color = INK2, linestyle = :dot, linewidth = 1)
        for (src, pk) in SOURCES
            G = grid_of(anchors(pk, ax), ax.logx)[1:2:end]
            g = [gmr_at(pk, ax, "Cmax", x) for x in G]
            sel = findall(x -> 0.70 <= x <= 1.40, g)
            y = [accept(pk, ax, "Joint", G[i], n; n_trial = N_PW) for i in sel]
            o = sortperm(g[sel])
            lines!(a, g[sel][o], y[o]; color = COL[src], linewidth = 2.5, label = src)
        end
    end
    Legend(fig[2, 1:2], [LineElement(color = COL[s], linewidth = 3) for s in first.(SOURCES)],
           first.(SOURCES); orientation = :horizontal, framevisible = false)
    save(joinpath(OUT, "fig_oc_cmax.png"), fig, px_per_unit = 2)
end

# (d) Type I error: PBPK truth at its own boundary vs at PopPK's boundary
let d = filter(r -> r.judged_by == "PBPK" && r.endpoint == "Joint" && !ismissing(r.T1_at_other_bnd), t1)
    if nrow(d) > 0
        labs = ["$(r.axis) $(r.side)\nn=$(r.n)" for r in eachrow(d)]
        fig = Figure(size = (1100, 480))
        a = Axis(fig[1, 1]; title = "Type I error (joint ABE), judged by the PBPK truth",
                 subtitle = "at the PBPK boundary (nominal) vs at the boundary PopPK would declare; dashed = 5%",
                 titlealign = :left, ylabel = "P(conclude BE)", xticks = (1:nrow(d), labs))
        barplot!(a, (1:nrow(d)) .- 0.19, d.T1_own; width = 0.36, color = COL["PBPK"],
                 label = "at PBPK boundary")
        barplot!(a, (1:nrow(d)) .+ 0.19, Float64.(d.T1_at_other_bnd); width = 0.36,
                 color = COL["PopPK"], label = "at PopPK boundary")
        hlines!(a, [0.05]; color = INK, linestyle = :dash)
        axislegend(a; position = :lt, framevisible = false)
        save(joinpath(OUT, "fig_type1_bars.png"), fig, px_per_unit = 2)
    end
end
println("\nWrote CSVs and figures to ", OUT)
