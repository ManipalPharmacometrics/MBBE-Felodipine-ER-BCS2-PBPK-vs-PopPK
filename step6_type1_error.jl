# =============================================================================
# STEP 6 — type I error at the BE boundary on each design axis, PBPK vs PopPK.
#
# In : results/step1_pbpk_endpoints.csv, results/step3_poppk_endpoints.csv
# Out: results/step6_boundaries.csv          where the true GMR hits 0.80 / 1.25
#      results/step6_type1_error.csv         pass rate at the boundary (target 5%)
#      results/step6_tost_control.csv        pass rate with the reference x 0.80
#      results/step6_scaled_at_boundary.csv  RSABE / ABEL pass rate at the ABE Cmax boundary
#      results/fig_step6_type1_bars.png, fig_step6_oc.png
#
# 1. Boundary: walk out from the reference along the axis (161-point grid)
#    until the noise-free GMR leaves 0.80–1.25 (Cmax and AUC separately;
#    Joint = the one reached first).
# 2. Type I error = share of 20,000 simulated BE trials that pass AT the boundary.
# 3. Cross-check: the pass rate at the boundary the OTHER method found.
#    PBPK is the truth, so "PBPK at the PopPK boundary" is the type I error a
#    PopPK-based safe space would really have.
# 4. Control: test = reference x 0.80 exactly (checks the TOST itself).
# =============================================================================

using CSV, DataFrames, CairoMakie
include("settings.jl")
include("be_functions.jl")

pbpk  = CSV.read(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), DataFrame)
poppk = CSV.read(joinpath(RESULTS, "step3_poppk_endpoints.csv"), DataFrame)
SOURCES = ["PBPK" => pbpk, "PopPK" => poppk]
COLOR = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

function wilson(p, n; z = 1.96)
    d = 1 + z^2 / n; c = (p + z^2 / (2n)) / d
    h = z * sqrt(p * (1 - p) / n + z^2 / (4n^2)) / d
    (max(0.0, c - h), min(1.0, c + h))
end

# ---- 1. boundaries ------------------------------------------------------------
bnd = DataFrame()
for ax in AXES, (source, pk) in SOURCES
    G = grid_of(pk, ax)
    gmr = Dict(ep => [gmr_at(pk, ax, ep, x) for x in G] for ep in ENDPOINTS)
    for side in (:low, :high)
        hits = Dict(ep => first_crossing(G, gmr[ep], ax.ref, side) for ep in ENDPOINTS)
        for ep in ENDPOINTS
            push!(bnd, (axis = ax.name, source = source, endpoint = ep, side = String(side),
                        boundary = hits[ep].x, limit = hits[ep].limit, binding = ep); promote = true)
        end
        reached = [ep for ep in ENDPOINTS if !ismissing(hits[ep].x)]
        first_ep = isempty(reached) ? nothing : reached[argmin([abs(hits[ep].x - ax.ref) for ep in reached])]
        push!(bnd, (axis = ax.name, source = source, endpoint = "Joint", side = String(side),
                    boundary = first_ep === nothing ? missing : hits[first_ep].x,
                    limit = first_ep === nothing ? "not reached" : hits[first_ep].limit,
                    binding = first_ep === nothing ? "none" : first_ep); promote = true)
    end
end
CSV.write(joinpath(RESULTS, "step6_boundaries.csv"), bnd)

# ---- 2–3. type I error at own and at the other method's boundary ----------------
t1 = DataFrame()
for ax in AXES, ep in [ENDPOINTS; "Joint"], side in ("low", "high"), n in N_SUBJ, (source, pk) in SOURCES
    sel(s) = (bnd.axis .== ax.name) .& (bnd.endpoint .== ep) .& (bnd.side .== side) .& s
    own   = only(bnd[sel(bnd.source .== source), :])
    other = only(bnd[sel(bnd.source .!= source), :])
    ismissing(own.boundary) && continue
    p = pass_rate_at(pk, ax, ep, own.boundary, n, N_TRIAL_TYPE1)
    lo, hi = wilson(p, N_TRIAL_TYPE1)
    p_other = ismissing(other.boundary) ? missing : pass_rate_at(pk, ax, ep, other.boundary, n, N_TRIAL_TYPE1)
    g_other = ismissing(other.boundary) ? missing : gmr_at(pk, ax, own.binding, other.boundary)
    push!(t1, (axis = ax.name, endpoint = ep, side = side, n_subj = n, judged_by = source,
               own_boundary = own.boundary, limit = own.limit, binding = own.binding,
               type1_own = p, ci_low = lo, ci_high = hi, other_boundary = other.boundary,
               type1_at_other_boundary = p_other, true_GMR_at_other_boundary = g_other); promote = true)
end
CSV.write(joinpath(RESULTS, "step6_type1_error.csv"), t1)

# ---- 4. TOST control (reference x 0.80) -------------------------------------------
ctl = DataFrame()
for (source, pk) in SOURCES, ep in [ENDPOINTS; "Joint"], n in N_SUBJ
    rC, rA = subject_values(pk, REF.T50, REF.shape, "Cmax"), subject_values(pk, REF.T50, REF.shape, "AUC_tEnd")
    p = ep == "Joint" ?
        pass_rate_joint(rC .* BE_LOWER, rC, rA .* BE_LOWER, rA; n_subj = n, n_trial = N_TRIAL_TYPE1,
                        cvC = WSV_CV["Cmax"], cvA = WSV_CV["AUC_tEnd"]) :
        pass_rate(subject_values(pk, REF.T50, REF.shape, ep) .* BE_LOWER, subject_values(pk, REF.T50, REF.shape, ep);
                  n_subj = n, n_trial = N_TRIAL_TYPE1, cv = WSV_CV[ep])
    push!(ctl, (source = source, endpoint = ep, n_subj = n, pass_rate = p))
end
CSV.write(joinpath(RESULTS, "step6_tost_control.csv"), ctl)

# ---- scaled methods at the ABE Cmax boundary (n = 24) ------------------------------
sc = DataFrame()
for ax in AXES, (source, pk) in SOURCES, side in ("low", "high")
    b = only(bnd[(bnd.axis .== ax.name) .& (bnd.source .== source) .& (bnd.endpoint .== "Cmax") .& (bnd.side .== side), :])
    ismissing(b.boundary) && continue
    for method in (:fda_rsabe, :ema_abel)
        p = pass_rate_scaled(values_at(pk, ax, "Cmax", b.boundary), values_at(pk, ax, "Cmax", ax.ref);
                             n_subj = N_SCALED, n_trial = N_TRIAL_TYPE1, method)
        push!(sc, (axis = ax.name, source = source, side = side, boundary = b.boundary,
                   method = String(method), n_subj = N_SCALED, pass_rate = p))
    end
end
CSV.write(joinpath(RESULTS, "step6_scaled_at_boundary.csv"), sc)

# ---- figures ----------------------------------------------------------------------
d = t1[(t1.judged_by .== "PBPK") .& (t1.endpoint .== "Joint") .& .!ismissing.(t1.type1_at_other_boundary), :]
fig = Figure(size = (900, 450))
a = Axis(fig[1, 1]; title = "Type I error, joint ABE, judged by PBPK", ylabel = "P(conclude BE)",
         xticks = (1:nrow(d), ["$(r.axis) $(r.side)\nn=$(r.n_subj)" for r in eachrow(d)]))
barplot!(a, (1:nrow(d)) .- 0.2, d.type1_own; width = 0.4, color = COLOR["PBPK"], label = "at PBPK boundary")
barplot!(a, (1:nrow(d)) .+ 0.2, Float64.(d.type1_at_other_boundary); width = 0.4,
         color = COLOR["PopPK"], label = "at PopPK boundary")
hlines!(a, [0.05]; color = :black, linestyle = :dash)
axislegend(a; position = :lt)
save(joinpath(RESULTS, "fig_step6_type1_bars.png"), fig)

ax = AXES[1]    # operating characteristic on the T50 axis (the shape axis never reaches 0.80)
fig = Figure(size = (1100, 450))
for (j, n) in enumerate(N_SUBJ)
    a = Axis(fig[1, j]; title = "Joint ABE, n = $n, T50 axis", xlabel = "true Cmax GMR", ylabel = "P(conclude BE)")
    vlines!(a, [BE_LOWER, BE_UPPER]; color = :gray60, linestyle = :dash)
    hlines!(a, [0.05]; color = :gray60, linestyle = :dot)
    for (source, pk) in SOURCES
        G = grid_of(pk, ax)[1:2:end]
        g = [gmr_at(pk, ax, "Cmax", x) for x in G]
        keep = findall(x -> 0.70 <= x <= 1.40, g)
        p = [pass_rate_at(pk, ax, "Joint", G[i], n, N_TRIAL_POWER) for i in keep]
        o = sortperm(g[keep])
        lines!(a, g[keep][o], p[o]; color = COLOR[source], label = source)
    end
    j == 1 && axislegend(a; position = :rt)
end
save(joinpath(RESULTS, "fig_step6_oc.png"), fig)

println("Boundaries:"); show(bnd; allrows = true); println()
println("\nType I error:")
show(t1[:, [:axis, :endpoint, :side, :n_subj, :judged_by, :own_boundary, :type1_own,
            :other_boundary, :type1_at_other_boundary]]; allrows = true); println()
