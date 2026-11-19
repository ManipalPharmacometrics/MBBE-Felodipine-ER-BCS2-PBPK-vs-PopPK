# =============================================================================
# STEP 5 — power for every formulation and along the two design axes,
# PBPK vs PopPK.
#
# In : results/step1_pbpk_endpoints.csv, results/step3_poppk_endpoints.csv
# Out: results/step5_power_by_formulation.csv  per formulation:
#        ABE 2x2, n = 54: pass rate for Cmax, AUC_tEnd and joint (both)
#        FDA RSABE and EMA ABEL, replicate design, n = 24, Cmax
#      results/step5_power_axes.csv            joint ABE power along each axis, n = 24 / 54
#      results/fig_step5_power_axes.png
#
# Pass rate = share of 1,000 simulated BE trials that pass (within-subject
# CVs in settings.jl). A formulation is inside the safe space if power >= 80%.
# =============================================================================

using CSV, DataFrames, CairoMakie
include("settings.jl")
include("be_functions.jl")

pbpk  = CSV.read(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), DataFrame)
poppk = CSV.read(joinpath(RESULTS, "step3_poppk_endpoints.csv"), DataFrame)
SOURCES = ["PBPK" => pbpk, "PopPK" => poppk]
COLOR = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

# ---- per formulation ----------------------------------------------------------
byform = DataFrame()
for (source, pk) in SOURCES, f in eachrow(sort(unique(pk[:, [:T50, :shape]]), [:shape, :T50]))
    test(ep) = subject_values(pk, f.T50, f.shape, ep)
    ref(ep)  = subject_values(pk, REF.T50, REF.shape, ep)
    abe(ep)  = pass_rate(test(ep), ref(ep); n_subj = N_ABE, n_trial = N_TRIAL_POWER, cv = WSV_CV[ep])
    push!(byform, (source = source, T50 = f.T50, shape = f.shape,
                   ABE_Cmax = abe("Cmax"), ABE_AUC = abe("AUC_tEnd"),
                   ABE_joint = pass_rate_joint(test("Cmax"), ref("Cmax"), test("AUC_tEnd"), ref("AUC_tEnd");
                                               n_subj = N_ABE, n_trial = N_TRIAL_POWER,
                                               cvC = WSV_CV["Cmax"], cvA = WSV_CV["AUC_tEnd"]),
                   RSABE_Cmax = pass_rate_scaled(test("Cmax"), ref("Cmax"); n_subj = N_SCALED,
                                                 n_trial = N_TRIAL_POWER, method = :fda_rsabe),
                   ABEL_Cmax  = pass_rate_scaled(test("Cmax"), ref("Cmax"); n_subj = N_SCALED,
                                                 n_trial = N_TRIAL_POWER, method = :ema_abel)))
end
byform.in_safe_space = byform.ABE_joint .>= POWER_TARGET
CSV.write(joinpath(RESULTS, "step5_power_by_formulation.csv"), byform)

# ---- along the axes -----------------------------------------------------------
axes_power = DataFrame()
for ax in AXES, (source, pk) in SOURCES, n in N_SUBJ, x in grid_of(pk, ax)[1:2:end]
    push!(axes_power, (axis = ax.name, source = source, n_subj = n, x = x,
                       power = pass_rate_at(pk, ax, "Joint", x, n, N_TRIAL_POWER)))
end
CSV.write(joinpath(RESULTS, "step5_power_axes.csv"), axes_power)

fig = Figure(size = (1100, 750))
for (i, n) in enumerate(reverse(N_SUBJ)), (j, ax) in enumerate(AXES)
    a = Axis(fig[i, j]; title = "Joint ABE power, n = $n", xlabel = ax.name, ylabel = "power",
             xscale = ax.log_x ? log10 : identity)
    hlines!(a, [POWER_TARGET]; color = :gray60, linestyle = :dash)
    vlines!(a, [ax.ref]; color = :gray60, linestyle = :dot)
    for source in ("PBPK", "PopPK")
        d = axes_power[(axes_power.axis .== ax.name) .& (axes_power.source .== source) .& (axes_power.n_subj .== n), :]
        lines!(a, d.x, d.power; color = COLOR[source], label = source)
    end
    i == 1 && j == 1 && axislegend(a; position = :rt)
end
save(joinpath(RESULTS, "fig_step5_power_axes.png"), fig)

println("Power by formulation (ABE n = $N_ABE; RSABE/ABEL n = $N_SCALED):")
show(byform; allrows = true, allcols = true); println()
