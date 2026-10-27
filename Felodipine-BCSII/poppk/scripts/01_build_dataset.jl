# =============================================================================
# 01_build_dataset.jl — felodipine ER: PBPK per-subject endpoints + simulated
# rich-sampling BE dataset for the PopPK fits. Mirrors the metformin 01.
#
# Truth = the 27 PK-Sim population runs (100 paired subjects, 0-96 h) in
#   Felodipine-BCSII/pksim-releaserate-cqa-er/matrix_Sweep_td50_shape/
# (T50 x shape cross + 2-D fill; reference T50 4.78 h, shape 0.764).
#
# Outputs (output/):
#   pbpk_pk_parameters_individual.csv  PK-Sim NCA per subject: T50, shape,
#                                      IndividualId, Cmax (nM), AUC_tEnd (nM·h, 0-96 h)
#   be_dataset.csv                     NONMEM-style dataset, 18 samples/subject,
#                                      15% proportional residual error
# Design: rich single-dose ER schedule to 96 h (sampling to >= 3 half-lives);
# residual error, seed and reference choice as in the metformin baseline.
# =============================================================================

using CSV, DataFrames, Random, Statistics

const PKDIR = normpath(joinpath(@__DIR__, "../../pksim-releaserate-cqa-er/matrix_Sweep_td50_shape"))
include(normpath(joinpath(@__DIR__, "../../../programs/julia/be_utils.jl"))); using .BEUtils

const HERE = @__DIR__
const OUT  = get(ENV, "POPPK_OUT", normpath(joinpath(HERE, "..", "output"))); mkpath(OUT)

const MW        = 384.26                    # felodipine g/mol
const DOSE_NMOL = 10.0 / MW * 1e6           # 10 mg -> 26 024 nmol
const RICH  = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 16, 24, 36, 48, 72, 96]
const SPARSE = [1.0, 3, 6, 12, 24, 48]
const SCHEDULE = get(ENV, "POPPK_SCHEDULE", "rich") == "sparse" ? SPARSE : Float64.(RICH)
const PROP_CV = 0.15
const SEED    = 20261016
const CONC_COL = "Organism|PeripheralVenousBlood|Felodipine|Plasma (Peripheral Venous Blood) [µmol/l]"

# ---- PBPK endpoints (PK-Sim NCA) --------------------------------------------
pk = load_all_formulations(PKDIR; drug = :felodipine)
pb = select(pk, :T50, :shape, :IndividualId,
            :Cmax => (x -> x .* 1e3) => :Cmax,                 # µmol/L -> nM
            :AUC_tEnd => (x -> x .* 1e3 ./ 60) => :AUC_tEnd)  # µmol·min/L -> nM·h
sort!(pb, [:T50, :shape, :IndividualId])
CSV.write(joinpath(OUT, "pbpk_pk_parameters_individual.csv"), pb)

# ---- simulated BE study from the dense profiles -----------------------------
files = sort(filter(f -> endswith(f, "-Results.csv") && startswith(f, "Pop_"), readdir(PKDIR)))
rng  = MersenneTwister(SEED)
rows = DataFrame(SUBJID = String[], T50 = Float64[], shape = Float64[], IndividualId = Int[],
                 time = Float64[], evid = Int[], amt = Union{Missing,Float64}[],
                 cmt = Int[], dv = Union{Missing,Float64}[])
for f in files
    fm = BEUtils.parse_formulation(f, :felodipine)
    df = CSV.read(joinpath(PKDIR, f), DataFrame)
    tm = df[!, "Time [min]"]; cc = df[!, CONC_COL]; ii = df[!, "IndividualId"]
    for id in sort(unique(ii))
        sel = ii .== id; st = tm[sel]; sc = cc[sel]
        sid = "T$(fm.T50)_b$(fm.shape)_ID$(id)"
        push!(rows, (sid, fm.T50, fm.shape, id, 0.0, 1, DOSE_NMOL, 1, missing))
        for th in SCHEDULE
            k = findfirst(==(round(th * 60)), st)
            k === nothing && error("t=$th h missing in $f")
            c = Float64(sc[k]) * 1e3 * (1 + PROP_CV * randn(rng))
            push!(rows, (sid, fm.T50, fm.shape, id, th, 0, 0.0, 2, max(c, 1e-6)))
        end
    end
end
sort!(rows, [:T50, :shape, :IndividualId, :time])
CSV.write(joinpath(OUT, "be_dataset.csv"), rows)
println("Formulations: ", nrow(unique(rows[:, [:T50, :shape]])), "  PBPK rows: ", nrow(pb))
println("Dose $(round(DOSE_NMOL, digits=0)) nmol; schedule ($(length(SCHEDULE)) obs): ", SCHEDULE)
println("Wrote be_dataset.csv ($(nrow(rows)) rows) and pbpk_pk_parameters_individual.csv -> $OUT")
