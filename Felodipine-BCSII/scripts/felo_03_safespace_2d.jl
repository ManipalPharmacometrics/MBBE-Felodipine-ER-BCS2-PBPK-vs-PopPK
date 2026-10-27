# =============================================================================
# felo_03_safespace_2d.jl
# Felodipine ER (MBBE drug two): 2-D T50 x shape safe-space figure.
#
# Input : <sweep dir> written by felo_02_be_sweep.jl (be_point_estimates.csv,
#         vbe_abe_2x2.csv, vbe_rsabe_cmax.csv).
# Method: power is a steep function of the GMR, so it is NOT interpolated
#         directly (that splits the band into islands on this sparse design).
#         Instead:
#           1. log GMR (Cmax, AUC_0-t; 100 paired subjects) is interpolated
#              linearly on a Delaunay triangulation of the simulated points, in
#              (log T50, shape) scaled coordinates. Nothing outside the hull.
#           2. log GMR -> power through a monotone (isotonic) calibration curve
#              of simulated pass rate vs |log GMR| over all formulations.
#           3. Joint ABE power = P_Cmax x P_AUC (checked against the simulated
#              joint rate at every point; max difference printed).
#         RSABE / ABEL contours use the Cmax criterion only (as in the sweep).
# Output: <sweep dir>/fig_safespace_2d_joint.png   filled joint ABE power,
#                     0.80 contours for ABE (n=54), RSABE and ABEL (n=24)
#         <sweep dir>/safespace_2d_boundaries.csv  T50 range with power >= 0.80
#                     at each simulated shape, per method; open_* = bound sits
#                     on the edge of the design, the true one lies beyond
#
# Overlay mode (2026-10-16): a second sweep dir (e.g. the PopPK endpoints run
# through felo_02) adds its joint-ABE 0.80 contour in orange and writes
# <sweep dir>/fig_safespace_2d_overlay_<label>.png instead of the single figure.
# AUC is AUC_0-t (as in the sweep's joint); until 2026-10-16 this used AUC_inf,
# which moved no boundary by more than the grid step (AUC is never binding).
#
# Usage (from programs, for the Julia env):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_03_safespace_2d.jl <sweep dir> [<overlay sweep dir> <label>]
# =============================================================================

using CSV, DataFrames, Printf, CairoMakie
# DelaunayTriangulation ships with Makie but is not a direct dep of this env.
const DT = Base.require(Base.PkgId(
    Base.UUID("927a84f5-c5f4-47a5-9785-b46e178433df"), "DelaunayTriangulation"))

const SWEEP   = ARGS[1]
const REF     = (4.78, 0.764)
const TARGET  = 0.80
const XLIM    = (3.3, 7.3)        # T50, h
const YLIM    = (0.55, 1.25)      # Weibull shape; above 1.2 only one point (1.5) supports it

mean(x) = sum(x) / length(x)

"Power surfaces for one sweep dir (felo_02 outputs). See the header for the method."
function build_surface(SWEEP)
    pe  = CSV.read(joinpath(SWEEP, "be_point_estimates.csv"), DataFrame)
    abe = CSV.read(joinpath(SWEEP, "vbe_abe_2x2.csv"), DataFrame)
    rs  = CSV.read(joinpath(SWEEP, "vbe_rsabe_cmax.csv"), DataFrame)
    n_abe   = only(unique(abe.n_subj))
    n_rsabe = only(unique(rs.n_subj))

    gmr(m) = select(filter(:metric => ==(m), pe), :T50, :shape, :GMR => (g -> log.(g)) => :lg)
    d = innerjoin(rename(gmr("Cmax"), :lg => :lgC), rename(gmr("AUC_tEnd"), :lg => :lgA),
                  on = [:T50, :shape])
    d = innerjoin(d, select(abe, :T50, :shape, :pass_Cmax, :pass_AUC_tEnd, :pass_joint),
                  on = [:T50, :shape])
    for (m, col) in (("fda_rsabe", :p_rsabe), ("ema_abel", :p_abel))
        d = innerjoin(d, select(filter(:method => ==(m), rs), :T50, :shape, :pass => col),
                             on = [:T50, :shape])
    end

    # --- 1. linear interpolation of log GMR on the triangulation -------------------
    scale(t, b) = (log(t / REF[1]) / 0.2, (b - REF[2]) / 0.3)
    pts  = [scale(t, b) for (t, b) in zip(d.T50, d.shape)]
    tris = [Tuple(T) for T in DT.each_solid_triangle(DT.triangulate(pts))]

    function interp(z, t, b)
        x, y = scale(t, b)
        for (i, j, k) in tris
            (x1, y1), (x2, y2), (x3, y3) = pts[i], pts[j], pts[k]
            den = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
            l1 = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / den
            l2 = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / den
            l3 = 1 - l1 - l2
            min(l1, l2, l3) >= -1e-9 && return l1 * z[i] + l2 * z[j] + l3 * z[k]
        end
        return NaN
    end

    # --- 2. monotone calibration |log GMR| -> pass rate ----------------------------
    function isotonic_decreasing(x, y)
        o = sortperm(x); xs = x[o]; v = Float64.(y[o]); w = ones(length(v))
        blocks = [[i] for i in eachindex(v)]
        i = 1
        while i < length(v)
            if v[i] < v[i + 1]                       # violation of "decreasing"
                m = (v[i] * w[i] + v[i + 1] * w[i + 1]) / (w[i] + w[i + 1])
                v[i] = m; w[i] += w[i + 1]; append!(blocks[i], blocks[i + 1])
                deleteat!(v, i + 1); deleteat!(w, i + 1); deleteat!(blocks, i + 1)
                i = max(i - 1, 1)
            else
                i += 1
            end
        end
        fx = [mean(xs[bk]) for bk in blocks]
        a -> a <= fx[1] ? v[1] : a >= fx[end] ? v[end] :
             (k = searchsortedlast(fx, a); v[k] + (v[k+1] - v[k]) * (a - fx[k]) / (fx[k+1] - fx[k]))
    end

    calC     = isotonic_decreasing(abs.(d.lgC), d.pass_Cmax)
    calA     = isotonic_decreasing(abs.(d.lgA), d.pass_AUC_tEnd)
    calRSABE = isotonic_decreasing(abs.(d.lgC), d.p_rsabe)
    calABEL  = isotonic_decreasing(abs.(d.lgC), d.p_abel)

    # --- 3. power surfaces ---------------------------------------------------------
    jk, rk, ek = "ABE (2x2, n=$n_abe)", "FDA RSABE (n=$n_rsabe)", "EMA ABEL (n=$n_rsabe)"
    function power(k, t, b)
        lc = interp(d.lgC, t, b); isnan(lc) && return NaN
        k == jk && return calC(abs(lc)) * calA(abs(interp(d.lgA, t, b)))
        k == rk && return calRSABE(abs(lc))
        return calABEL(abs(lc))
    end
    chk = [abs(power(jk, r.T50, r.shape) - r.pass_joint) for r in eachrow(d)]
    @printf("Calibrated vs simulated joint ABE power at the %d points: max |diff| %.3f, median %.3f\n",
            nrow(d), maximum(chk), sort(chk)[cld(length(chk), 2)])

    xs = range(XLIM..., length = 241)
    ys = range(YLIM..., length = 201)
    grids = Dict(k => [power(k, x, y) for x in xs, y in ys] for k in (jk, rk, ek))
    (d = d, abe = abe, power = power, grids = grids, xs = xs, ys = ys,
     n_abe = n_abe, n_rsabe = n_rsabe, jk = jk, rk = rk, ek = ek)
end

S = build_surface(SWEEP)
OV = length(ARGS) >= 3 ? (dir = ARGS[2], label = ARGS[3], s = build_surface(ARGS[2])) : nothing
d, abe, power, grids, xs, ys = S.d, S.abe, S.power, S.grids, S.xs, S.ys
n_abe, n_rsabe, jk, rk, ek = S.n_abe, S.n_rsabe, S.jk, S.rk, S.ek

# --- boundaries at each simulated shape --------------------------------------
bounds = DataFrame(method = String[], shape = Float64[], T50_lo = Float64[],
                   T50_hi = Float64[], ratio_lo = Float64[], ratio_hi = Float64[],
                   open_lo = Bool[], open_hi = Bool[])
for k in (jk, rk, ek), b in sort(unique(d.shape))
    fi = (t, b) -> power(k, t, b)
    ts = range(2.4, 9.5, step = 0.005)
    ok = [t for t in ts if (v = fi(t, b); !isnan(v) && v >= TARGET)]
    isempty(ok) && continue
    # open_* = the bound sits on the edge of the design, so the true one lies beyond
    lo, hi = first(ok), last(ok)
    push!(bounds, (k, b, lo, hi, lo / REF[1], hi / REF[1],
                   isnan(fi(lo - 0.01, b)), isnan(fi(hi + 0.01, b))))
end
sort!(bounds, [:method, :shape])
OV === nothing && CSV.write(joinpath(SWEEP, "safespace_2d_boundaries.csv"), bounds)
println(bounds)

# --- figure -------------------------------------------------------------------
ink, ink2 = colorant"#1f1f1d", colorant"#5f5e5a"
blues = [colorant"#cde2fb", colorant"#86b6ef", colorant"#3987e5", colorant"#184f95"]
levels = [0.0, 0.5, 0.8, 0.9, 1.0]

fig = Figure(size = (900, 680), fontsize = 15, backgroundcolor = :white)
Label(fig[0, 1:2], "Felodipine ER 10 mg: where is a test formulation bioequivalent?";
      halign = :left, font = :bold, fontsize = 18)
Label(fig[1, 1:2], (OV === nothing ? "" : "Fill and black lines: PBPK (PK-Sim). Orange: $(OV.label). ") *
      "Fill: joint ABE power (Cmax and AUC₀₋ₜ, 2x2, n=$n_abe). Lines: power = 0.80 " *
      "per method (RSABE/ABEL: Cmax, n=$n_rsabe). Between simulated points: interpolated GMR, " *
      "calibrated to power."; halign = :left, justification = :left, color = ink2, fontsize = 13,
      word_wrap = true)
ax = Axis(fig[2, 1];
    xlabel = "Dissolution T50 (h)", ylabel = "Weibull shape",
    xgridcolor = (:black, 0.06), ygridcolor = (:black, 0.06),
    limits = (XLIM, YLIM))
top = Axis(fig[2, 1]; xaxisposition = :top, xlabel = "T50 / reference T50",
    limits = (XLIM ./ REF[1], YLIM), xgridvisible = false, ygridvisible = false,
    yticksvisible = false, yticklabelsvisible = false, xlabelcolor = ink2,
    xticklabelcolor = ink2)
hidespines!(top); top.xticks = 0.7:0.1:1.5

cf = contourf!(ax, xs, ys, grids[jk]; levels = levels, colormap = blues)
for (k, (ls, lw)) in [jk => (:solid, 3), rk => (:dash, 2), ek => (:dot, 2)]
    contour!(ax, xs, ys, grids[k]; levels = [TARGET], color = ink,
             linestyle = ls, linewidth = lw, label = "$k = 0.80")
end

if OV !== nothing
    contour!(ax, xs, ys, OV.s.grids[OV.s.jk]; levels = [TARGET], color = colorant"#eb6834",
             linewidth = 3.5)
    # contour! legend entries ignore `color`; a NaN line carries the legend entry
    lines!(ax, [NaN], [NaN]; color = colorant"#eb6834", linewidth = 3.5,
           label = "$(OV.label): ABE (2x2, n=$(OV.s.n_abe)) = 0.80")
end

pass = abe.pass_joint .>= TARGET
inview = (XLIM[1] .<= abe.T50 .<= XLIM[2]) .& (YLIM[1] .<= abe.shape .<= YLIM[2])
scatter!(ax, abe.T50[pass .& inview], abe.shape[pass .& inview];
         color = ink, markersize = 12, strokewidth = 0,
         label = "simulated, joint power ≥ 0.80")
scatter!(ax, abe.T50[.!pass .& inview], abe.shape[.!pass .& inview];
         color = :white, markersize = 11, strokecolor = ink, strokewidth = 1.5,
         label = "simulated, joint power < 0.80")
scatter!(ax, [REF[1]], [REF[2]]; marker = :star5, markersize = 20,
         color = OV === nothing ? colorant"#eb6834" : :white,
         strokecolor = OV === nothing ? :white : ink, strokewidth = 1.5,
         label = "reference (Plendil ER, OSP Weibull)")

Colorbar(fig[2, 2], cf; label = "joint ABE power", ticks = levels)
Legend(fig[3, 1:2], ax; orientation = :horizontal, nbanks = 2, framevisible = false,
       labelsize = 13)
fname = OV === nothing ? "fig_safespace_2d_joint.png" :
        "fig_safespace_2d_overlay_$(replace(OV.label, r"[^A-Za-z0-9]" => "")).png"
save(joinpath(SWEEP, fname), fig, px_per_unit = 2)
println("Wrote ", joinpath(SWEEP, fname))
