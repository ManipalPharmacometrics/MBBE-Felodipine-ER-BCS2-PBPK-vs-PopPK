# =============================================================================
# 03_collect.jl — gather the per-arm PopPK results from 02 into the tables the
# comparison uses, and copy everything to the learning project.
#
# Outputs (output/ and LEARN_DST):
#   poppk_pk_parameters_individual.csv  T50, shape, IndividualId, Cmax (nM),
#                                       AUC_tEnd (nM·h, 0-96 h), AUC_inf (nM·h), Tmax
#   poppk_param_estimates.csv, poppk_param_se.csv, poppk_shrinkage.csv
#   poppk_ipred_thin.csv                 0.5-h IPRED grid, all arms
#   pk_analyses_format/Pop_T50_<h>h_b<shape>-PK-Analyses.csv
#       the PopPK endpoints in PK-Sim's PK-Analyses layout and units, so the
#       PBPK sweep scripts (felo_02/felo_03) run on them unchanged.
# AUC_inf = AUC_tEnd + C(96 h)/lambda_z, lambda_z from the 72-96 h IPRED
# (log-linear; exact for the 1-cmt model's terminal phase).
# =============================================================================

using CSV, DataFrames, Statistics

const ROOT = normpath(joinpath(@__DIR__, ".."))
const OUT  = get(ENV, "POPPK_OUT", joinpath(ROOT, "output"))
const ARMS = joinpath(OUT, "arms")
const LEARN_DST = get(ENV, "LEARN_DST", normpath(joinpath(@__DIR__, "../../Felodipine-BCSII/poppk/output_rich")))

dirs = sort(filter(d -> isfile(joinpath(ARMS, d, "nca.csv")), readdir(ARMS)))
cat(name) = reduce(vcat, [CSV.read(joinpath(ARMS, d, name), DataFrame)
                          for d in dirs if isfile(joinpath(ARMS, d, name))]; cols = :union)
nca = cat("nca.csv"); thin = cat("ipred_thin.csv")
par = cat("params.csv"); shr = cat("shrinkage.csv")
se  = try; cat("se.csv"); catch; DataFrame(); end

lz = combine(groupby(thin, [:T50, :shape, :IndividualId])) do g
    s = g[g.time .>= 72, :]
    k = -(log(s.ipred[end]) - log(s.ipred[1])) / (s.time[end] - s.time[1])
    (lambda_z = k, Clast = s.ipred[end])
end
nca = leftjoin(nca, lz, on = [:T50, :shape, :IndividualId])
nca.AUC_inf = nca.AUC_tEnd .+ nca.Clast ./ nca.lambda_z
sort!(nca, [:T50, :shape, :IndividualId])

mkpath(LEARN_DST)
for dst in (OUT, LEARN_DST)
    CSV.write(joinpath(dst, "poppk_pk_parameters_individual.csv"),
              select(nca, :T50, :shape, :IndividualId, :Cmax, :AUC_tEnd, :AUC_inf, :Tmax))
    CSV.write(joinpath(dst, "poppk_param_estimates.csv"), sort(par, [:shape, :T50]))
    CSV.write(joinpath(dst, "poppk_shrinkage.csv"), sort(shr, [:shape, :T50]))
    nrow(se) > 0 && CSV.write(joinpath(dst, "poppk_param_se.csv"), se)
    CSV.write(joinpath(dst, "poppk_ipred_thin.csv"), thin)
    pkd = joinpath(dst, "pk_analyses_format"); mkpath(pkd)
    for g in groupby(nca, [:T50, :shape])
        rows = DataFrame(IndividualId = Int[], QuantityPath = String[], Parameter = String[],
                         Value = Float64[], Unit = String[])
        qp = "Organism|PeripheralVenousBlood|Felodipine|Plasma (Peripheral Venous Blood)"
        for r in eachrow(g)
            push!(rows, (r.IndividualId, qp, "C_max", r.Cmax / 1e3, "µmol/l"))
            push!(rows, (r.IndividualId, qp, "AUC_tEnd", r.AUC_tEnd * 60 / 1e3, "µmol*min/l"))
            push!(rows, (r.IndividualId, qp, "AUC_inf", r.AUC_inf * 60 / 1e3, "µmol*min/l"))
        end
        CSV.write(joinpath(pkd, "Pop_T50_$(g.T50[1])h_b$(g.shape[1])-PK-Analyses.csv"), rows)
    end
end
# PBPK table and the dataset travel with the results, for provenance
for f in ("pbpk_pk_parameters_individual.csv", "be_dataset.csv")
    cp(joinpath(OUT, f), joinpath(LEARN_DST, f); force = true)
end
println("Collected $(length(dirs)) arms -> $OUT and $LEARN_DST")
showcols = intersect([:T50, :shape, :basin, :T50_implied, :tvt50, :tvgamma, :tvka, :tvcl, :tvvc, :tvq, :tvvp,
                      :sigma_p, :LL, :LL_other], Symbol.(names(par)))   # v1 (1-cmt) or v2 columns
show(select(sort(par, [:shape, :T50]), showcols); allrows = true, allcols = true); println()
