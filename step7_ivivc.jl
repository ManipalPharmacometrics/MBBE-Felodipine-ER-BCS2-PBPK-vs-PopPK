# =============================================================================
# STEP 7 — IVIVC check: does the PopPK release follow the PK-Sim release?
#
# In : results/step2_parameters.csv   (M5 estimates for the 27 formulations)
# Out: results/step7_ivivc.csv        per formulation: PK-Sim T50/shape vs PopPK
#      results/step7_ivivc_fit.csv    regression slopes along the two design axes
#      results/fig_step7_ivivc.png
#
# How PopPK relates to PBPK here:
#   PK-Sim T50 and shape = the Weibull release set in each simulation (the truth).
#   M5 releases through TWO Weibull components, so it is summarised by the
#   release curve as a whole (lag excluded):
#       F(t) = FF * (1 - 2^-(t/T50a)^Ga) + (1 - FF) * (1 - 2^-(t/T50b)^Gb)
#     effective T50   = time at which F = 50%
#     effective shape = the single-Weibull shape through F = 25% and 75%:
#                       log(log 4 / log(4/3)) / log(t75 / t25)
#   T50 axis  (shape 0.764): slope of log(PopPK T50) on log(PK-Sim T50); IVIVC -> 1
#   shape axis (T50 4.78 h): slope of PopPK shape on PK-Sim shape;        IVIVC -> 1
# =============================================================================

using CSV, DataFrames, Statistics, Distributions, CairoMakie
include("settings.jl")

p = CSV.read(joinpath(RESULTS, "step2_parameters.csv"), DataFrame)

weibull(t, t50, g) = 1 - 2.0^(-(t / t50)^g)
released(r, t) = r.tvff * weibull(t, r.tvt50a, r.tvga) + (1 - r.tvff) * weibull(t, r.tvt50b, r.tvgb)
"Time at which the M5 release curve reaches fraction q (bisection)."
function time_to(r, q)
    lo, hi = 0.0, 1.0
    while released(r, hi) < q; hi *= 2; end
    for _ in 1:60; mid = (lo + hi) / 2; released(r, mid) < q ? (lo = mid) : (hi = mid); end
    (lo + hi) / 2
end

ivivc = DataFrame(T50_PKSim = p.T50, shape_PKSim = p.shape,
                  T50_PopPK = [time_to(r, 0.5) for r in eachrow(p)],
                  shape_PopPK = [log(log(4) / log(4 / 3)) / log(time_to(r, 0.75) / time_to(r, 0.25)) for r in eachrow(p)],
                  fraction_fast = p.tvff, CL = p.tvcl)
CSV.write(joinpath(RESULTS, "step7_ivivc.csv"), ivivc)

"Ordinary least-squares slope with a 95% CI."
function ols(x, y)
    b = cov(x, y) / var(x); a = mean(y) - b * mean(x); n = length(x)
    se = sqrt(sum((y .- a .- b .* x) .^ 2) / (n - 2) / sum((x .- mean(x)) .^ 2))
    tq = quantile(TDist(n - 2), 0.975)
    (slope = b, ci_lower = b - tq * se, ci_upper = b + tq * se)
end
on_T50   = ivivc[ivivc.shape_PKSim .== REF.shape, :]
on_shape = ivivc[ivivc.T50_PKSim .== REF.T50, :]
fits = DataFrame([
    merge((analysis = "T50 axis: log(PopPK T50) on log(PK-Sim T50)",), ols(log.(on_T50.T50_PKSim), log.(on_T50.T50_PopPK))),
    merge((analysis = "shape axis: PopPK shape on PK-Sim shape",), ols(on_shape.shape_PKSim, on_shape.shape_PopPK))])
fits.ivivc_would_give .= 1.0
CSV.write(joinpath(RESULTS, "step7_ivivc_fit.csv"), fits)

fig = Figure(size = (1100, 480))
a = Axis(fig[1, 1]; title = "T50 axis (shape 0.764)", xlabel = "PK-Sim T50 (h)", ylabel = "PopPK effective T50 (h)",
         xscale = log10, yscale = log10)
lines!(a, [2.3, 10.0], [2.3, 10.0]; color = :gray60, linestyle = :dash, label = "identity (perfect IVIVC)")
scatterlines!(a, sort(on_T50, :T50_PKSim).T50_PKSim, sort(on_T50, :T50_PKSim).T50_PopPK; color = :royalblue, label = "M5 fits")
axislegend(a; position = :lt)
a = Axis(fig[1, 2]; title = "Shape axis (T50 4.78 h)", xlabel = "PK-Sim shape", ylabel = "PopPK effective shape")
lines!(a, [0.45, 1.55], [0.45, 1.55]; color = :gray60, linestyle = :dash)
scatterlines!(a, sort(on_shape, :shape_PKSim).shape_PKSim, sort(on_shape, :shape_PKSim).shape_PopPK; color = :royalblue)
save(joinpath(RESULTS, "fig_step7_ivivc.png"), fig)

println("PopPK vs PK-Sim release:"); show(sort(ivivc, [:shape_PKSim, :T50_PKSim]); allrows = true); println()
println("\nSlopes (1 = PopPK release follows PK-Sim):"); show(fits; allcols = true); println()
