# =============================================================================
# felo_05_poppk_diagnostics.jl
# Felodipine ER PopPK (per-arm, 1-cmt Weibull release): how well the fitted
# model reproduces the PK-Sim truth it was fitted to.
#
#   fig_profiles.png    mean PK-Sim profile vs mean PopPK IPRED, 5 arms (0-24 h
#                       linear + 0-96 h log inset), with the simulated BE samples
#   fig_endpoint_agreement.png  per-subject PopPK / PBPK ratio of Cmax and AUC_0-t,
#                       all 27 arms (box per arm, ordered by T50 within shape)
#   fig_release_params.png  fitted release (implied T50, shape) vs the PK-Sim input
#   endpoint_agreement.csv  GM and 5-95% of the per-subject ratios per arm
#
# Inputs : <dir>/{pbpk_pk_parameters_individual, poppk_pk_parameters_individual,
#          poppk_param_estimates, poppk_ipred_thin, be_dataset}.csv
#          + the PK-Sim Results.csv files (dense truth)
# Usage (from programs):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_05_poppk_diagnostics.jl ../Felodipine-BCSII/poppk/output_rich
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie, Printf

const DIR   = ARGS[1]
const OUT   = joinpath(DIR, "diagnostics"); mkpath(OUT)
const PKDIR = joinpath(@__DIR__, "..", "pksim-releaserate-cqa-er", "matrix_Sweep_td50_shape")
const CONC  = "Organism|PeripheralVenousBlood|Felodipine|Plasma (Peripheral Venous Blood) [µmol/l]"
COL = Dict("PBPK" => colorant"#2a78d6", "PopPK" => colorant"#eb6834")
INK, INK2 = colorant"#1f1f1d", colorant"#5f5e5a"
set_theme!(fontsize = 15, Axis = (xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
                                  topspinevisible = false, rightspinevisible = false))
gm(v) = exp(mean(log.(v)))

pb  = CSV.read(joinpath(DIR, "pbpk_pk_parameters_individual.csv"), DataFrame)
pp  = CSV.read(joinpath(DIR, "poppk_pk_parameters_individual.csv"), DataFrame)
par = CSV.read(joinpath(DIR, "poppk_param_estimates.csv"), DataFrame)
ip  = CSV.read(joinpath(DIR, "poppk_ipred_thin.csv"), DataFrame)
obs = filter(r -> r.evid == 0, CSV.read(joinpath(DIR, "be_dataset.csv"), DataFrame))

# ---- per-subject endpoint agreement ---------------------------------------------
j = innerjoin(pb, pp, on = [:T50, :shape, :IndividualId], renamecols = "_pb" => "_pp")
j.rC = j.Cmax_pp ./ j.Cmax_pb; j.rA = j.AUC_tEnd_pp ./ j.AUC_tEnd_pb
agr = combine(groupby(j, [:T50, :shape]),
    :rC => gm => :Cmax_GM, :rC => (x -> quantile(x, 0.05)) => :Cmax_p05,
    :rC => (x -> quantile(x, 0.95)) => :Cmax_p95,
    :rA => gm => :AUC_GM, :rA => (x -> quantile(x, 0.05)) => :AUC_p05,
    :rA => (x -> quantile(x, 0.95)) => :AUC_p95)
sort!(agr, [:shape, :T50])
CSV.write(joinpath(OUT, "endpoint_agreement.csv"), agr)
@printf("Per-subject PopPK/PBPK: Cmax GM %.3f (arms %.3f-%.3f), AUC_0-t GM %.3f (arms %.3f-%.3f)\n",
        gm(j.rC), extrema(agr.Cmax_GM)..., gm(j.rA), extrema(agr.AUC_GM)...)

let fig = Figure(size = (1300, 720))
    Label(fig[0, 1], "PopPK individual predictions vs the PK-Sim truth, per subject (100 per arm)";
          font = :bold, fontsize = 17, halign = :left, tellwidth = false)
    ord = sort(unique(j[:, [:T50, :shape]]), [:shape, :T50])
    labs = ["$(r.T50) h · $(r.shape)" for r in eachrow(ord)]
    pos = Dict((r.T50, r.shape) => i for (i, r) in enumerate(eachrow(ord)))
    for (row, (col, lab)) in enumerate(((:rC, "Cmax"), (:rA, "AUC₀₋ₜ")))
        a = Axis(fig[row, 1]; ylabel = "PopPK / PBPK, $lab",
                 xticks = (1:nrow(ord), labs), xticklabelsize = 11, xticklabelrotation = π / 2,
                 xticklabelsvisible = row == 2,
                 xlabel = row == 2 ? "formulation: T50 · shape" : "")
        hlines!(a, [1.0]; color = INK2, linestyle = :dash)
        hspan!(a, 0.95, 1.05; color = (:black, 0.05))
        boxplot!(a, [pos[(r.T50, r.shape)] for r in eachrow(j)], j[!, col];
                 color = COL["PopPK"], show_outliers = false, whiskerwidth = 0.5, width = 0.6)
        ylims!(a, 0.8, 1.2)
    end
    save(joinpath(OUT, "fig_endpoint_agreement.png"), fig, px_per_unit = 2)
end

# ---- fitted release parameters vs the PK-Sim input ---------------------------
let fig = Figure(size = (1200, 480))
    Label(fig[0, 1:2], "Fitted PopPK release vs the PK-Sim formulation input (one fit per arm)";
          font = :bold, fontsize = 17, halign = :left)
    a1 = Axis(fig[1, 1]; xlabel = "PK-Sim dissolution T50 (h)", ylabel = "PopPK implied T50 (h)",
              title = "all arms", titlealign = :left)
    a2 = Axis(fig[1, 2]; xlabel = "PK-Sim Weibull shape", ylabel = "PopPK Weibull shape",
              title = "all arms", titlealign = :left)
    v2 = hasproperty(par, :tvt50)              # v2 (Weibull -> ka -> 2-cmt) estimates T50 directly
    for (a, x, y) in ((a1, par.T50, v2 ? par.tvt50 : par.T50_implied), (a2, par.shape, par.tvgamma))
        lo, hi = extrema(vcat(x, y))
        lines!(a, [lo, hi], [lo, hi]; color = INK2, linestyle = :dash)
        scatter!(a, x, y; color = [b == "release_limited" ? COL["PopPK"] : :white for b in par.basin],
                 strokecolor = COL["PopPK"], strokewidth = 1.5, markersize = 11)
    end
    Legend(fig[2, 1:2], [MarkerElement(marker = :circle, color = :white, strokecolor = COL["PopPK"], strokewidth = 1.5),
                         MarkerElement(marker = :circle, color = COL["PopPK"])],
           v2 ? ["absorption-limited start (higher likelihood)", "release-limited start (higher likelihood)"] :
                ["disposition-limited fit (higher likelihood)", "release-limited fit (higher likelihood)"];
           orientation = :horizontal, framevisible = false)
    save(joinpath(OUT, "fig_release_params.png"), fig, px_per_unit = 2)
end

# ---- mean profiles: PK-Sim truth vs PopPK IPRED vs simulated samples ------------
ARMS = [(2.39, 0.764), (4.78, 0.764), (9.56, 0.764), (4.78, 0.5), (4.78, 1.5)]
include(joinpath(@__DIR__, "..", "..", "programs", "julia", "be_utils.jl"))
const RESULTS = Dict(Tuple(BEUtils.parse_formulation(f, :felodipine)) => joinpath(PKDIR, f)
                     for f in readdir(PKDIR) if startswith(f, "Pop_") && endswith(f, "-Results.csv"))
let fig = Figure(size = (1500, 560))
    Label(fig[0, 1:length(ARMS)], "Mean plasma profile: PK-Sim truth, PopPK fit, and the simulated BE samples (n = 100 each)";
          font = :bold, fontsize = 17, halign = :left)
    for (k, (t, s)) in enumerate(ARMS)
        df = CSV.read(RESULTS[(t, s)], DataFrame)
        tr = combine(groupby(df, "Time [min]"), CONC => (c -> mean(c) * 1e3) => :c)
        tr.h = tr[!, "Time [min]"] ./ 60
        pr = combine(groupby(filter(r -> r.T50 == t && r.shape == s, ip), :time), :ipred => mean => :c)
        ob = combine(groupby(filter(r -> r.T50 == t && r.shape == s, obs), :time), :dv => mean => :c)
        for (row, (xmax, ysc)) in enumerate(((24.0, identity), (96.0, log10)))
            a = Axis(fig[row, k]; xlabel = row == 2 ? "time (h)" : "",
                     ylabel = k == 1 ? (row == 1 ? "felodipine (nM)" : "nM, log") : "",
                     title = row == 1 ? "T50 $(t) h, shape $(s)" : "", titlealign = :left,
                     yscale = ysc, titlesize = 14,
                     yticks = row == 2 ? ([0.1, 1, 10], ["0.1", "1", "10"]) : Makie.automatic)
            xlims!(a, 0, xmax)
            sel(d, col) = d[(d[!, col] .<= xmax) .& (d.c .> (row == 2 ? 0.05 : 0)), :]
            a1 = sel(tr, :h); a2 = sel(pr, :time); a3 = sel(ob, :time)
            lines!(a, a1.h, a1.c; color = COL["PBPK"], linewidth = 3)
            lines!(a, a2.time, a2.c; color = COL["PopPK"], linewidth = 2, linestyle = :dash)
            scatter!(a, a3.time, a3.c; color = INK, markersize = 6)
            row == 1 ? ylims!(a, 0, 13) : ylims!(a, 0.05, 20)
        end
    end
    Legend(fig[3, 1:length(ARMS)], [LineElement(color = COL["PBPK"], linewidth = 3),
                                     LineElement(color = COL["PopPK"], linewidth = 2, linestyle = :dash),
                                     MarkerElement(marker = :circle, color = INK, markersize = 7)],
           ["PK-Sim truth (mean)", "PopPK IPRED (mean)", "simulated BE samples (mean, 15% error)"];
           orientation = :horizontal, framevisible = false)
    save(joinpath(OUT, "fig_profiles.png"), fig, px_per_unit = 2)
end
println("Wrote diagnostics to ", OUT)
