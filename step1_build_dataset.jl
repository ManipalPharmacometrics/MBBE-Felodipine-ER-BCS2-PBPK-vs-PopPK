# =============================================================================
# STEP 1 — PBPK endpoints and the PopPK dataset from the PK-Sim exports.
#
# In : data/pksim_exports/Pop_<formulation>-Results.csv      plasma profiles,
#                                                             100 subjects, 0–96 h
#      data/pksim_exports/Pop_<formulation>-PK-Analyses.csv  PK-Sim's own NCA
# Out: results/step1_pbpk_endpoints.csv   T50, shape, IndividualId, Cmax (nM),
#                                         AUC_tEnd (nM·h, 0–96 h) — the PBPK "truth"
#      results/step1_dataset.csv          NONMEM-style: one dose row + 18 samples
#                                         per subject per formulation
#
# Formulation from the file name: Pop_T50_<h>h_b<shape>; a missing token means
# the reference value (T50 4.78 h, shape 0.764), e.g. Pop_b0.5 = 4.78 h / 0.5.
# Sampled concentrations get 15% proportional assay error:
#     obs = true * (1 + 0.15 * N(0,1))
# =============================================================================

using CSV, DataFrames, Random
include("settings.jl")
mkpath(RESULTS)

const PLASMA = "Organism|PeripheralVenousBlood|Felodipine|Plasma (Peripheral Venous Blood)"

function formulation(file)
    t = match(r"T50_([0-9.]+)h", file)
    b = match(r"_b([0-9]*\.[0-9]+)", file)
    (T50 = t === nothing ? REF.T50 : parse(Float64, t.captures[1]),
     shape = b === nothing ? REF.shape : parse(Float64, b.captures[1]))
end

files(kind) = sort(filter(f -> endswith(f, "-$kind.csv"), readdir(DATA)))

# ---- 1. PBPK endpoints from PK-Sim's NCA (µmol/L -> nM, µmol·min/L -> nM·h) --
pbpk = DataFrame()
for file in files("PK-Analyses")
    f   = formulation(file)
    nca = CSV.read(joinpath(DATA, file), DataFrame)
    nca = nca[nca.QuantityPath .== PLASMA, :]
    cmax = sort(nca[nca.Parameter .== "C_max", :], :IndividualId)
    auc  = sort(nca[nca.Parameter .== "AUC_tEnd", :], :IndividualId)
    append!(pbpk, DataFrame(T50 = f.T50, shape = f.shape, IndividualId = cmax.IndividualId,
                            Cmax = cmax.Value .* 1e3, AUC_tEnd = auc.Value .* 1e3 ./ 60))
end
sort!(pbpk, [:T50, :shape, :IndividualId])
CSV.write(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), pbpk)

# ---- 2. sampled dataset for the PopPK fits ----------------------------------
rng  = MersenneTwister(SEED_DATASET)
rows = DataFrame()
for file in files("Results")                                  # file-name order = noise order
    f   = formulation(file)
    sim = CSV.read(joinpath(DATA, file), DataFrame)
    for id in sort(unique(sim.IndividualId))
        subj = sim[sim.IndividualId .== id, :]
        sid  = "T$(f.T50)_b$(f.shape)_ID$(id)"
        push!(rows, (SUBJID = sid, T50 = f.T50, shape = f.shape, IndividualId = id, time = 0.0,
                     evid = 1, amt = DOSE_NMOL, cmt = 1, dv = missing); promote = true)
        for t in SCHEDULE
            k = findfirst(==(round(t * 60)), subj[!, "Time [min]"])      # PK-Sim time is in minutes
            c = subj[k, PLASMA * " [µmol/l]"] * 1e3 * (1 + PROP_CV * randn(rng))   # nM + error
            push!(rows, (SUBJID = sid, T50 = f.T50, shape = f.shape, IndividualId = id, time = float(t),
                         evid = 0, amt = 0.0, cmt = 2, dv = max(c, 1e-6)); promote = true)
        end
    end
end
sort!(rows, [:T50, :shape, :IndividualId, :time])
CSV.write(joinpath(RESULTS, "step1_dataset.csv"), rows)

println("Formulations: ", nrow(unique(pbpk[:, [:T50, :shape]])), "   subjects each: ",
        length(unique(pbpk.IndividualId)))
println("Samples per subject: ", length(SCHEDULE), "  (", join(SCHEDULE, ", "), " h)")
println("Wrote results/step1_pbpk_endpoints.csv and results/step1_dataset.csv (", nrow(rows), " rows)")
