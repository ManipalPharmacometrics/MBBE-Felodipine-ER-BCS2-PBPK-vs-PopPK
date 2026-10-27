# =============================================================================
# felo_08_curve_tmax_check.jl — does the PopPK curve shape (absorption phase,
# Tmax) match PK-Sim, by model variant? (2026-10-24)
#
# Compares, on the diagnostic arms, the best start of each fitted model in
# poppk-fits/output_modelcomp/ with the noise-free
# PK-Sim profiles (3-min grid). Written to answer whether adding the release lag
# to the @delay (@delay-form) models improves the curve and Tmax:
#   M1, M2, M3        Weibull hazard ODE, lag 0.3113 h
#   M1d, M2d, M3d     @delay form, lag 0.3113 h
#   M1d0, M2d0, M3d0  @delay form, no lag
#   M3L, M5, M6       structural variants for the early Tmax (06_structure_fits.jl):
#                     M3 with estimated lag, dual Weibull, 3-transit chain
# Extra fit folders can be added with FELO08_EXTRA="label=dir;label=dir"
# (dir holds <model>/<arm>/<start>/ like output_modelcomp).
#
# Outputs (poppk/model_comparison/curve_tmax/):
#   tmax_by_model.csv     per model x arm: median Tmax (model vs PK-Sim), median
#                         per-subject Tmax difference, % of subjects within 0.5 h,
#                         mean-curve error at 0.5/1/2/3/4/6 h and over 0.5-12 h
#   fig_mean_curves.png   mean IPRED vs mean PK-Sim, 0-12 h, one panel per arm
#                         (IPRED is stored on a 0.5-h grid, so 0-0.5 h is a straight line)
# Usage (from programs):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_08_curve_tmax_check.jl
# =============================================================================

using CSV, DataFrames, Statistics, CairoMakie, Printf
include(joinpath(@__DIR__, "..", "..", "programs", "julia", "be_utils.jl")); using .BEUtils

const FITS  = normpath(joinpath(@__DIR__, "../../poppk-fits"))
const PKDIR = joinpath(@__DIR__, "..", "pksim-releaserate-cqa-er", "matrix_Sweep_td50_shape")
const OUT   = joinpath(@__DIR__, "..", "poppk", "model_comparison", "curve_tmax"); mkpath(OUT)
const CONC  = "Organism|PeripheralVenousBlood|Felodipine|Plasma (Peripheral Venous Blood) [µmol/l]"
const EARLY = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0]
const ORDER = ["M1", "M1d", "M1d0", "M2", "M2d", "M2d0", "M3", "M3d", "M3d0", "M3L", "M5", "M6"]

# ---- fitted models: best start per model x arm --------------------------------
sources = [("", joinpath(FITS, "output_modelcomp"))]
for kv in split(get(ENV, "FELO08_EXTRA", ""), ";", keepempty = false)
    l, d = split(kv, "="); push!(sources, (String(l), String(d)))
end
best = Dict{Tuple{String,String},String}()          # (model, arm) -> start dir
for (lab, root) in sources, m in filter(x -> isdir(joinpath(root, x)), readdir(root)), a in readdir(joinpath(root, m))
    (m in ORDER || lab != "") || continue
    m == "M4" && continue
    ll = -Inf; bd = ""
    for st in readdir(joinpath(root, m, a))
        f = joinpath(root, m, a, st, "fitstats.csv"); isfile(f) || continue
        v = CSV.read(f, DataFrame).LL[1]
        (v isa Number && v > ll) && (ll = v; bd = joinpath(root, m, a, st))
    end
    bd == "" || (best[(m * lab, a)] = bd)
end
arms = sort(unique(last.(collect(keys(best)))))

# ---- PK-Sim truth -------------------------------------------------------------
pk = Dict{String,DataFrame}()                        # arm -> IndividualId, time (h), conc (nM)
for f in filter(f -> startswith(f, "Pop_") && endswith(f, "-Results.csv"), readdir(PKDIR))
    fm = BEUtils.parse_formulation(f, :felodipine); a = "$(fm.T50)_$(fm.shape)"
    a in arms || continue
    d = CSV.read(joinpath(PKDIR, f), DataFrame)
    pk[a] = DataFrame(IndividualId = d.IndividualId, time = d[!, "Time [min]"] ./ 60, conc = d[!, CONC] .* 1e3)
end
tmax_of(t, c) = t[argmax(c)]
interp(t, c, x) = (k = searchsortedlast(t, x); k >= length(t) ? c[end] :
                   c[k] + (c[k+1] - c[k]) * (x - t[k]) / (t[k+1] - t[k]))

meancurve(df, grid) = [mean(interp(sdf.time, sdf.conc, x) for sdf in groupby(sort(df, [:IndividualId, :time]), :IndividualId))
                       for x in grid]

rows = DataFrame(); curves = Dict{Tuple{String,String},Vector{Float64}}()
const GRID = collect(0.0:0.25:12.0)
for a in arms
    haskey(pk, a) || continue
    p = pk[a]
    ptm = combine(groupby(p, :IndividualId), [:time, :conc] => ((t, c) -> tmax_of(t, c)) => :Tmax_pk)
    curves[("PK-Sim", a)] = meancurve(p, GRID)
    for m in sort([k[1] for k in keys(best) if k[2] == a], by = x -> something(findfirst(==(x), ORDER), 99))
        dir = best[(m, a)]
        ep = CSV.read(joinpath(dir, "endpoints.csv"), DataFrame; types = Dict(:period => String))
        th = CSV.read(joinpath(dir, "ipred_thin.csv"), DataFrame; types = Dict(:period => String))
        th = rename(th, :ipred => :conc)
        j = innerjoin(ep[:, [:IndividualId, :Tmax]], ptm, on = :IndividualId)
        d = j.Tmax .- j.Tmax_pk
        mc = meancurve(th, GRID); curves[(m, a)] = mc
        pkc = curves[("PK-Sim", a)]
        early = Dict(Symbol(replace("meanC_err_pct_$(x)h", "." => "p")) => 100 * (mc[findfirst(==(x), GRID)] / pkc[findfirst(==(x), GRID)] - 1) for x in EARLY)
        sel = GRID .>= 0.5                        # PK-Sim is 0 before its 0.31 h lag
        push!(rows, merge((model = m, arm = a, start = basename(dir),
                           Tmax_median_model = median(j.Tmax), Tmax_median_pksim = median(j.Tmax_pk),
                           Tmax_meancurve_model = GRID[argmax(mc)], Tmax_meancurve_pksim = GRID[argmax(pkc)],
                           dTmax_median = median(d), dTmax_abs_median = median(abs.(d)),
                           pct_within_0p5h = 100 * mean(abs.(d) .<= 0.5),
                           meanC_abs_err_pct_0p5_12h = 100 * mean(abs.(mc[sel] ./ pkc[sel] .- 1))),
                          (; early...)); cols = :union)
    end
end
CSV.write(joinpath(OUT, "tmax_by_model.csv"), rows)

# ---- figure ---------------------------------------------------------------------
COL = Dict("PK-Sim" => :black, "M1" => colorant"#eb6834", "M1d" => colorant"#f2a07b", "M1d0" => colorant"#8a3a17",
           "M2" => colorant"#3aa35b", "M2d" => colorant"#8fd1a3", "M2d0" => colorant"#1e5c32",
           "M3" => colorant"#2a78d6", "M3d" => colorant"#8fb8ec", "M3d0" => colorant"#173f73",
           "M3L" => colorant"#7b4fc9", "M5" => colorant"#d6336c", "M6" => colorant"#c99a06")
fig = Figure(size = (1600, 420 * cld(length(arms), 3)), fontsize = 13)
for (i, a) in enumerate(arms)
    haskey(pk, a) || continue
    ax = Axis(fig[cld(i, 3), mod1(i, 3)], title = "arm $(replace(a, "_" => " h / shape "))",
              xlabel = "time (h)", ylabel = "mean conc (nM)")
    lines!(ax, GRID, curves[("PK-Sim", a)], color = :black, linewidth = 3, label = "PK-Sim")
    extra = sort([k[1] for k in keys(curves) if k[2] == a && !(k[1] in ORDER) && k[1] != "PK-Sim"])
    for m in vcat(ORDER, extra)
        haskey(curves, (m, a)) || continue
        nolag = endswith(m, "d0") || occursin("nolag", m)
        lines!(ax, GRID, curves[(m, a)], color = get(COL, m, get(COL, m[1:2] * "d0", :gray)), linewidth = 1.8,
               linestyle = nolag ? :dash : endswith(m, "d") ? :dot : :solid, label = m)
    end
    i == 1 && axislegend(ax, position = :rt, nbanks = 2)
end
save(joinpath(OUT, "fig_mean_curves.png"), fig, px_per_unit = 1.5)

show(select(rows, :arm, :model, :Tmax_median_model, :Tmax_median_pksim, :dTmax_median, :pct_within_0p5h,
            :meanC_err_pct_0p5h, :meanC_err_pct_1p0h, :meanC_err_pct_2p0h, :meanC_abs_err_pct_0p5_12h), allrows = true, allcols = true)
println("\nWrote -> $OUT")
