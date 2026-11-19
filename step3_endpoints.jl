# =============================================================================
# STEP 3 — BE endpoints (Cmax, AUC 0–96 h) per subject from the PopPK
# individual predictions, and a check against PK-Sim.
#
# In : results/fits/*/ipred.csv, results/step1_pbpk_endpoints.csv
# Out: results/step3_poppk_endpoints.csv   T50, shape, IndividualId, Cmax (nM), AUC_tEnd (nM·h)
#      results/step3_check_vs_pbpk.csv     geometric means per formulation, PopPK/PBPK ratio
#
#   Cmax     = NCA.cmax: highest predicted concentration (0.05 h grid)
#   AUC_tEnd = NCA.auc (auctype :last, method :linear): AUC 0–96 h, linear trapezoid
# =============================================================================

using CSV, DataFrames, Statistics, NCA
include("settings.jl")

geomean(x) = exp(mean(log.(x)))

poppk = DataFrame()
for dir in readdir(joinpath(RESULTS, "fits"); join = true)
    ipred = CSV.read(joinpath(dir, "ipred.csv"), DataFrame)
    # Pumas NCA package: one NCA subject per IndividualId (one formulation per folder)
    pop = read_nca(ipred; id = :IndividualId, time = :time, observations = :ipred,
                   group = [:T50, :shape], verbose = false)
    nca = innerjoin(NCA.cmax(pop), NCA.auc(pop; auctype = :last, method = :linear),
                    on = [:id, :T50, :shape])
    append!(poppk, DataFrame(T50 = nca.T50, shape = nca.shape, IndividualId = parse.(Int, string.(nca.id)),
                             Cmax = nca.cmax, AUC_tEnd = nca.auc))
end
sort!(poppk, [:T50, :shape, :IndividualId])
CSV.write(joinpath(RESULTS, "step3_poppk_endpoints.csv"), poppk)

pbpk = CSV.read(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), DataFrame)
gm(df, src) = combine(groupby(df, [:T50, :shape]), :Cmax => geomean => Symbol(src, "_Cmax"),
                      :AUC_tEnd => geomean => Symbol(src, "_AUC"))
check = innerjoin(gm(pbpk, "PBPK"), gm(poppk, "PopPK"), on = [:T50, :shape])
check.ratio_Cmax = check.PopPK_Cmax ./ check.PBPK_Cmax
check.ratio_AUC  = check.PopPK_AUC ./ check.PBPK_AUC
sort!(check, [:shape, :T50])
CSV.write(joinpath(RESULTS, "step3_check_vs_pbpk.csv"), check)
println("Geometric means (Cmax nM, AUC 0–96 h nM·h) and PopPK/PBPK ratios:")
show(check; allrows = true, allcols = true); println()
