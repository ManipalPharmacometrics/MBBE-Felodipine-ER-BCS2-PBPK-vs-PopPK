# =============================================================================
# STEP 4 — test/reference GMR with 90% CI for every formulation, PBPK vs PopPK.
#
# In : results/step1_pbpk_endpoints.csv, results/step3_poppk_endpoints.csv
# Out: results/step4_gmr.csv       all 27 formulations, both sources
#      results/fig_step4_gmr.png   GMR along the T50 axis and the shape axis
#
# Each formulation vs the reference (T50 4.78 h, shape 0.764) in the same 100
# virtual subjects, no within-subject noise: the "true" GMR. The test is the
# standard Pumas one, `pumas_be` (Bioequivalence package): a 2x2 crossover
# ANOVA, 90% CI of the GMR. Odd-numbered subjects get sequence RT, even TR.
# The reference row (reference vs itself) is GMR 1 by definition and has no CI.
# =============================================================================

using CSV, DataFrames, CairoMakie, Bioequivalence
include("settings.jl")
include("be_functions.jl")

pbpk  = CSV.read(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), DataFrame)
poppk = CSV.read(joinpath(RESULTS, "step3_poppk_endpoints.csv"), DataFrame)
COLOR = Dict("PBPK" => :firebrick, "PopPK" => :royalblue)

gmr = DataFrame()
for (source, pk) in ("PBPK" => pbpk, "PopPK" => poppk), ep in ENDPOINTS
    ref = subject_values(pk, REF.T50, REF.shape, ep)
    for f in eachrow(sort(unique(pk[:, [:T50, :shape]]), [:shape, :T50]))
        (f.T50, f.shape) == (REF.T50, REF.shape) && continue
        test = subject_values(pk, f.T50, f.shape, ep)
        n    = length(ref)
        seq  = [isodd(i) ? "RT" : "TR" for i in 1:n]
        crossover = DataFrame(id = repeat(1:n, inner = 2), sequence = repeat(seq, inner = 2),
                              period = repeat([1, 2], n))
        crossover.value = [s[p] == 'T' ? test[i] : ref[i]
                           for (i, s, p) in zip(crossover.id, crossover.sequence, crossover.period)]
        r = pumas_be(crossover; endpoint = :value, sigdigits = 16).result[1, :]
        push!(gmr, (source = source, endpoint = ep, T50 = f.T50, shape = f.shape, GMR = r.GMR,
                    CI_lower = r.LB, CI_upper = r.UB, pass = r.LB >= BE_LOWER && r.UB <= BE_UPPER))
    end
end
CSV.write(joinpath(RESULTS, "step4_gmr.csv"), gmr)

fig = Figure(size = (1100, 750))
for (i, ep) in enumerate(ENDPOINTS), (j, ax) in enumerate(AXES)
    a = Axis(fig[i, j]; title = "$ep GMR along the $(ax.name) axis", xlabel = ax.name,
             ylabel = "GMR vs reference", xscale = ax.log_x ? log10 : identity)
    hlines!(a, [BE_LOWER, BE_UPPER]; color = :gray60, linestyle = :dash)
    for source in ("PBPK", "PopPK")
        d = gmr[(gmr.source .== source) .& (gmr.endpoint .== ep) .& (gmr[!, ax.fixed] .== ax.fixed_value), :]
        sort!(d, ax.col)
        scatterlines!(a, d[!, ax.col], d.GMR; color = COLOR[source], label = source)
        rangebars!(a, d[!, ax.col], d.CI_lower, d.CI_upper; color = COLOR[source])
    end
    i == 1 && j == 1 && axislegend(a; position = :rt)
end
save(joinpath(RESULTS, "fig_step4_gmr.png"), fig)

show(unstack(gmr[:, [:source, :endpoint, :T50, :shape, :GMR]], [:endpoint, :T50, :shape], :source, :GMR);
     allrows = true); println()
