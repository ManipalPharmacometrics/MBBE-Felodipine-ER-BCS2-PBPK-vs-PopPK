# =============================================================================
# felo_02_be_sweep.jl
# Felodipine ER (MBBE drug two): BE across the coarse T50 x shape cross.
#
# Input : PK-Sim population exports (*-PK-Analyses.csv), 100 paired subjects,
#         reference = shipped OSP ER Weibull (T50 4.78 h, shape 0.764) with
#         colonic CYP3A4 at 50% of SI (project notes §9.6, §9.14-9.15).
# Output: <out>/be_point_estimates.csv   paired GMR + 90% CI, all 100 subjects
#         <out>/vbe_abe_2x2.csv          ABE pass rate, 2x2, n from power calc
#         <out>/vbe_rsabe_cmax.csv       Cmax reference-scaled sensitivity
#
# Usage (from programs, for the Julia env):
#   julia --project=. ../Felodipine-BCSII/scripts/felo_02_be_sweep.jl <pk dir> <out dir> [arm]
# arm = fasted (default) | fed. "fed" is a statistical sensitivity only: the
# fasted PK-Sim GMRs with the fed within-subject CVs and n. PK-Sim cannot
# reproduce the felodipine food effect (project notes §9.22), so there is no
# fed PK to sweep.
# =============================================================================

using CSV, DataFrames, Statistics, Printf
include(joinpath(@__DIR__, "..", "..", "programs", "julia", "be_utils.jl"))
using .BEUtils

const DRUG    = :felodipine
const PKDIR   = ARGS[1]
const OUT     = ARGS[2]
const ARM     = length(ARGS) >= 3 ? Symbol(ARGS[3]) : :fasted
const N_TRIAL = 1000
const N_RSABE = 24          # FDA minimum for a replicate RSABE study
const SEED    = 1234        # VBE resampling seed (be_utils default, stated explicitly)
mkpath(OUT)

pk   = load_all_formulations(PKDIR; drug = DRUG)
ref  = sort(reference_rows(pk, DRUG), :IndividualId)
rf   = ref_formulation(DRUG)
cv   = wsv_cv(DRUG, ARM)
cvwr = wsv_cv(DRUG, :cvwr_rsabe, "Cmax")
n2x2 = sample_size_abe_2x2(cv["Cmax"])          # 54 fasted (CV 0.357), 56 fed; GMR 0.95
forms = sort(unique(pk[:, [:T50, :shape]]), [:shape, :T50])
println("Formulations: ", nrow(forms), "   ref ", rf, "   arm ", ARM, "   n(2x2) = ", n2x2)

arm(t, s) = t == rf.T50 && s == rf.shape ? "ref" : s == rf.shape ? "T50" : "shape"

# ---- paired point estimates --------------------------------------------------
point = DataFrame()
for f in eachrow(forms)
    tst = sort(filter(r -> r.T50 == f.T50 && r.shape == f.shape, pk), :IndividualId)
    res = be_all_metrics(tst, ref)
    res[!, :T50] .= f.T50; res[!, :shape] .= f.shape; res[!, :axis] .= arm(f.T50, f.shape)
    append!(point, res)
end
CSV.write(joinpath(OUT, "be_point_estimates.csv"), point)

# ---- VBE: 2x2 ABE at the powered n ------------------------------------------
vbe = DataFrame()
for f in eachrow(forms)
    tst = sort(filter(r -> r.T50 == f.T50 && r.shape == f.shape, pk), :IndividualId)
    row = Dict{Symbol,Any}(:T50 => f.T50, :shape => f.shape, :axis => arm(f.T50, f.shape),
                           :n_subj => n2x2)
    for m in ("Cmax", "AUC_inf", "AUC_tEnd")
        r = vbe_replicate(tst[!, m], ref[!, m]; seed = SEED, n_subj = n2x2, n_trial = N_TRIAL, wsv_cv = cv[m])
        row[Symbol("pass_" * m)] = r.pass_rate
    end
    row[:pass_joint] = vbe_replicate_joint(tst.Cmax, ref.Cmax, tst.AUC_tEnd, ref.AUC_tEnd;
                                           seed = SEED, n_subj = n2x2, n_trial = N_TRIAL,
                                           cv_C = cv["Cmax"], cv_A = cv["AUC_tEnd"]).pass_rate
    push!(vbe, row; cols = :union)
end
select!(vbe, :axis, :T50, :shape, :n_subj, :pass_Cmax, :pass_AUC_tEnd, :pass_AUC_inf, :pass_joint)
CSV.write(joinpath(OUT, "vbe_abe_2x2.csv"), vbe)

# ---- Cmax reference-scaled sensitivity (full replicate) ----------------------
rs = DataFrame()
for f in eachrow(forms), method in (:fda_rsabe, :ema_abel)
    tst = sort(filter(r -> r.T50 == f.T50 && r.shape == f.shape, pk), :IndividualId)
    r = vbe_replicate_scaled(tst.Cmax, ref.Cmax; seed = SEED, n_subj = N_RSABE, n_trial = N_TRIAL,
                             cv_wr = cvwr, method)
    push!(rs, (axis = arm(f.T50, f.shape), T50 = f.T50, shape = f.shape, method = String(method),
               n_subj = N_RSABE, pass = r.pass_rate, scaled = r.scaled_rate))
end
CSV.write(joinpath(OUT, "vbe_rsabe_cmax.csv"), rs)

# ---- console summary -----------------------------------------------------------
cm = filter(:metric => ==("Cmax"), point); au = filter(:metric => ==("AUC_tEnd"), point)
fda = filter(:method => ==("fda_rsabe"), rs); ema = filter(:method => ==("ema_abel"), rs)
@printf("\n%-6s %5s %6s | %7s %7s | %6s %6s %6s | %6s %6s\n", "axis", "T50", "shape",
        "GMR_C", "GMR_A", "ABE_C", "ABE_A", "joint", "RSABE", "ABEL")
for (i, f) in enumerate(eachrow(forms))
    @printf("%-6s %5.2f %6.3f | %7.3f %7.3f | %6.3f %6.3f %6.3f | %6.3f %6.3f\n",
            arm(f.T50, f.shape), f.T50, f.shape, cm.GMR[i], au.GMR[i],
            vbe.pass_Cmax[i], vbe.pass_AUC_tEnd[i], vbe.pass_joint[i], fda.pass[i], ema.pass[i])
end
println("\nWrote CSVs to ", OUT)
