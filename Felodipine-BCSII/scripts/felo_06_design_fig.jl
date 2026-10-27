# felo_06_design_fig.jl — the 27-formulation PK-Sim design, by batch (deck figure).
# Usage (from programs): julia --project=. ../Felodipine-BCSII/scripts/felo_06_design_fig.jl
using CairoMakie
const OUT = joinpath(@__DIR__, "..", "figures_discussion"); mkpath(OUT)
coarse = [(t, 0.764) for t in (2.39, 3.35, 4.07, 5.74, 7.17, 9.56)] ∪
         [(4.78, b) for b in (0.5, 0.6, 0.9, 1.0, 1.2, 1.5)]
fill2d = [(t, b) for t in (4.07, 5.74) for b in (0.6, 0.9, 1.0, 1.2)] ∪ [(4.42, 0.764), (5.26, 0.764)]
extra  = [(6.31, 0.6), (6.31, 1.0), (6.31, 1.2), (6.88, 1.2)]
INK = colorant"#1f1f1d"
fig = Figure(size = (900, 560), fontsize = 15)
ax = Axis(fig[1, 1]; xlabel = "Dissolution T50 (h)", ylabel = "Weibull shape", xscale = log10,
          xticks = ([2.39, 3.35, 4.78, 7.17, 9.56], ["2.39", "3.35", "4.78", "7.17", "9.56"]),
          title = "27 PK-Sim population runs (100 paired subjects each, 0–96 h)", titlealign = :left,
          xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
          topspinevisible = false, rightspinevisible = false)
for (pts, c, m, lab) in ((coarse, colorant"#2a78d6", :circle, "coarse cross, 15 Oct (12)"),
                         (fill2d, colorant"#eb6834", :rect, "2-D fill, 16 Oct (10)"),
                         (extra,  colorant"#1baf7a", :diamond, "band extension, 16 Oct (4)"))
    scatter!(ax, first.(pts), last.(pts); color = c, marker = m, markersize = 15,
             strokecolor = :white, strokewidth = 1.5, label = lab)
end
scatter!(ax, [4.78], [0.764]; marker = :star5, markersize = 26, color = :white,
         strokecolor = INK, strokewidth = 2, label = "reference: Plendil ER (OSP Weibull)")
Legend(fig[2, 1], ax; orientation = :horizontal, nbanks = 2, framevisible = false)
save(joinpath(OUT, "fig_design.png"), fig, px_per_unit = 2)
println("Wrote ", joinpath(OUT, "fig_design.png"))
