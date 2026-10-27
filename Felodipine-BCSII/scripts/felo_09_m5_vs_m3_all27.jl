# =============================================================================
# felo_09_m5_vs_m3_all27.jl — final model M5 (dual Weibull @delay + lag) vs M3
# (fallback) on all 27 arms (2026-10-25).
#   per arm: LL/AIC (M3 12 est. par incl. no lag -> 11; M5 14 -> 14), Cmax and
#   AUC PopPK/PBPK GM ratio, T/R Cmax GMR vs the 4.78/0.764 reference (model vs
#   PBPK), median per-subject Tmax difference vs noise-free PK-Sim
#   and the 2-D safe-space boundaries (PBPK, M3, M5).
# Outputs: poppk/model_comparison/m5_vs_m3_all27/{per_arm.csv, boundaries.csv, summary.txt}
# Usage (from programs):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_09_m5_vs_m3_all27.jl
# =============================================================================

using CSV, DataFrames, Statistics, Printf
include(joinpath(@__DIR__, "..", "..", "programs", "julia", "be_utils.jl")); using .BEUtils

const P     = joinpath(@__DIR__, "..", "poppk")
const OUT   = joinpath(P, "model_comparison", "m5_vs_m3_all27"); mkpath(OUT)
const PKDIR = joinpath(@__DIR__, "..", "pksim-releaserate-cqa-er", "matrix_Sweep_td50_shape")
const CONC  = "Organism|PeripheralVenousBlood|Felodipine|Plasma (Peripheral Venous Blood) [µmol/l]"
const REF   = (4.78, 0.764)
const SRC   = Dict("M3" => joinpath(P, "output_rich_v2"), "M5" => joinpath(P, "output_rich_m5"))
const NPAR  = Dict("M3" => 11, "M5" => 14)

pb = CSV.read(joinpath(SRC["M5"], "pbpk_pk_parameters_individual.csv"), DataFrame)
gm(x) = exp(mean(log.(x)))

# PK-Sim noise-free Tmax per subject
pkt = DataFrame()
for f in filter(f -> startswith(f, "Pop_") && endswith(f, "-Results.csv"), readdir(PKDIR))
    fm = BEUtils.parse_formulation(f, :felodipine)
    d = CSV.read(joinpath(PKDIR, f), DataFrame)
    for g in groupby(d, :IndividualId)
        push!(pkt, (T50 = fm.T50, shape = fm.shape, IndividualId = g.IndividualId[1],
                    Tmax_pk = g[argmax(g[!, CONC]), "Time [min]"] / 60))
    end
end

rows = DataFrame()
for m in ("M3", "M5")
    par = CSV.read(joinpath(SRC[m], "poppk_param_estimates.csv"), DataFrame)
    pp  = CSV.read(joinpath(SRC[m], "poppk_pk_parameters_individual.csv"), DataFrame)
    j   = innerjoin(pp, select(pb, :T50, :shape, :IndividualId, :Cmax => :Cmax_pb, :AUC_tEnd => :AUC_pb),
                    on = [:T50, :shape, :IndividualId])
    j   = innerjoin(j, pkt, on = [:T50, :shape, :IndividualId])
    ref = sort(filter(r -> (r.T50, r.shape) == REF, j), :IndividualId)
    for g in groupby(j, [:T50, :shape])
        g = sort(g, :IndividualId); pr = par[(par.T50 .== g.T50[1]) .& (par.shape .== g.shape[1]), :]
        gmr_m = gm(g.Cmax ./ ref.Cmax); gmr_p = gm(g.Cmax_pb ./ ref.Cmax_pb)
        push!(rows, (model = m, T50 = g.T50[1], shape = g.shape[1], LL = pr.LL[1],
                     AIC = -2pr.LL[1] + 2NPAR[m],
                     Cmax_ratio = gm(g.Cmax ./ g.Cmax_pb), AUC_ratio = gm(g.AUC_tEnd ./ g.AUC_pb),
                     GMR_Cmax_model = gmr_m, GMR_Cmax_PBPK = gmr_p, GMR_err_pct = 100 * (gmr_m / gmr_p - 1),
                     dTmax_median = median(g.Tmax .- g.Tmax_pk), Tmax_pksim_median = median(g.Tmax_pk),
                     pct_Tmax_within_0p5h = 100 * mean(abs.(g.Tmax .- g.Tmax_pk) .<= 0.5)))
    end
end
sort!(rows, [:shape, :T50, :model]); CSV.write(joinpath(OUT, "per_arm.csv"), rows)

w = innerjoin(select(rows[rows.model .== "M3", :], :T50, :shape, :AIC => :AIC_M3, :GMR_err_pct => :err_M3,
                     :dTmax_median => :dT_M3, :pct_Tmax_within_0p5h => :w_M3),
              select(rows[rows.model .== "M5", :], :T50, :shape, :AIC => :AIC_M5, :GMR_err_pct => :err_M5,
                     :dTmax_median => :dT_M5, :pct_Tmax_within_0p5h => :w_M5), on = [:T50, :shape])
w.dAIC = w.AIC_M5 .- w.AIC_M3

b = DataFrame()
for (lab, d) in (("PBPK", joinpath(@__DIR__, "..", "output_sweep_2d")), ("M3", joinpath(P, "output_sweep_2d_poppk_v2")),
                 ("M5", joinpath(P, "output_sweep_2d_poppk_m5")))
    x = CSV.read(joinpath(d, "safespace_2d_boundaries.csv"), DataFrame); x.source .= lab; append!(b, x; cols = :union)
end
CSV.write(joinpath(OUT, "boundaries.csv"), b)

open(joinpath(OUT, "summary.txt"), "w") do io
    for o in (stdout, io)
        show(o, w, allrows = true, allcols = true); println(o)
        for m in ("M3", "M5")
            r = rows[rows.model .== m, :]
            @printf(o, "%s: |GMR err| max %.1f%% mean %.1f%%; Cmax ratio %.3f-%.3f; AUC ratio %.3f-%.3f; |dTmax| median-of-arms %.2f h, max %.2f h; within 0.5 h mean %.0f%%\n",
                    m, maximum(abs.(r.GMR_err_pct)), mean(abs.(r.GMR_err_pct)), extrema(r.Cmax_ratio)..., extrema(r.AUC_ratio)...,
                    median(abs.(r.dTmax_median)), maximum(abs.(r.dTmax_median)), mean(r.pct_Tmax_within_0p5h))
        end
        @printf(o, "M5 better AIC on %d/27 arms; dAIC median %.1f, range %.1f to %.1f\n", sum(w.dAIC .< 0), median(w.dAIC), extrema(w.dAIC)...)
        show(o, sort(select(b, :source, :method, :shape, :T50_lo, :T50_hi), [:method, :shape, :source]), allrows = true); println(o)
    end
end
