# =============================================================================
# STEP 9 — one tested formulation: predict BE for all the others (felodipine).
#
# In : results/step1_dataset.csv, results/step1_pbpk_endpoints.csv,
#      results/step3_poppk_endpoints.csv, model_known_release.jl
# Out: results/one_test/step9_parameters.csv, step9_predicted_endpoints.csv,
#      step9_gmr.csv, step9_decisions.csv, step9_boundaries.csv, fig_step9_gmr.png
#
# Scenario: a BE study has been run with the reference (T50 4.78 h, shape
# 0.764) and ONE test formulation. Can a PopPK model fitted to those two arms
# predict the other 25 formulations from their dissolution curves alone?
#   test "T50"  : T50 7.17 h, shape 0.764 (slower release, on the T50 axis)
#   test "shape": T50 4.78 h, shape 1.2   (different shape, on the shape axis)
# Model: model_known_release.jl — release = the known PK-Sim dissolution curve
# (not estimated); estimated = absorption, regional gut availability and
# 2-compartment disposition.
# tau (time scale of the move from small intestine to colon) is FIXED at 2, 3.5
# and 5 h, spanning typical small-intestinal transit: estimated freely it ran
# to ~0.2 h, where availability and CL/volumes cannot be separated.
# Method names: "T50_tau3.5" = tested on the T50 axis, tau 3.5 h.
# Prediction: each subject keeps the random effects of its reference arm; only
# the dissolution curve changes. Same 100 subjects as PK-Sim (paired BE).
# Compared with PBPK (truth) and with the M5 fits of step 2, which used data
# from all 27 formulations (the best PopPK can do when every arm is tested).
# =============================================================================

using Pumas, CSV, DataFrames, CairoMakie, NCA
include("settings.jl")
include("be_functions.jl")
include("model_known_release.jl")

const TESTED = ["T50" => (T50 = 7.17, shape = 0.764), "shape" => (T50 = 4.78, shape = 1.2)]
const OUT    = joinpath(RESULTS, "one_test"); mkpath(OUT)
const PRED_TIMES = collect(0.0:0.05:96.0)
const TAU_FIXED  = [2.0, 3.5, 5.0]               # h

data = CSV.read(joinpath(RESULTS, "step1_dataset.csv"), DataFrame)
data.cmt .= 1                                    # the dose goes into Tablet
pbpk = CSV.read(joinpath(RESULTS, "step1_pbpk_endpoints.csv"), DataFrame)
m5   = CSV.read(joinpath(RESULTS, "step3_poppk_endpoints.csv"), DataFrame)
FORMS = sort(unique(data[:, [:T50, :shape]]), [:shape, :T50])
to_population(df) = read_pumas(df; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                               observations = [:dv], covariates = [:T50, :shape, :IndividualId])
is_form(df, f) = (df.T50 .== f.T50) .& (df.shape .== f.shape)

# ---- 1. fit on the reference + one tested formulation, predict all 27 ----------------
refrows = data[is_form(data, REF), :]
params = DataFrame(); predicted = Dict{String,DataFrame}()
for (axis_name, tested) in TESTED, tau in TAU_FIXED
    name  = "$(axis_name)_tau$(tau)"
    train = to_population(data[is_form(data, REF) .| is_form(data, tested), :])
    @info "fitting $name (tested $(tested.T50) h, shape $(tested.shape))"
    fitted = fit(poppk_known_release, train, merge(init_params(poppk_known_release), (tvtau = tau,)), FOCE();
                 constantcoef = (:tvtau,))
    est = coef(fitted)
    push!(params, merge((method = name, tested_T50 = tested.T50, tested_shape = tested.shape,
                         minus2LL = -2 * loglikelihood(fitted)),
                        (; (k => v for (k, v) in pairs(est) if v isa Number)...),
                        (omega_tau = est.Ω[1, 1], omega_cl = est.Ω[2, 2], omega_vc = est.Ω[3, 3])))

    ebe = Dict(subj.id => η for (subj, η) in zip(train, empirical_bayes(fitted)))
    prof = DataFrame()
    for f in eachrow(FORMS)
        rows = copy(refrows); rows.T50 .= f.T50; rows.shape .= f.shape     # same subject, new formulation
        pop  = to_population(rows)
        for (subj, pred) in zip(pop, predict(poppk_known_release, pop, est, [ebe[s.id] for s in pop];
                                             obstimes = PRED_TIMES))
            append!(prof, DataFrame(T50 = f.T50, shape = f.shape,
                                    IndividualId = Int(subj.covariates(0.0).IndividualId),
                                    time = pred.time, ipred = pred.ipred.dv))
        end
    end
    # Cmax and AUC 0–96 h with Pumas NCA, as in step 3
    pop_nca = read_nca(prof; id = :IndividualId, time = :time, observations = :ipred,
                       group = [:T50, :shape], verbose = false)
    nca = innerjoin(NCA.cmax(pop_nca), NCA.auc(pop_nca; auctype = :last, method = :linear),
                    on = [:id, :T50, :shape])
    predicted[name] = DataFrame(T50 = nca.T50, shape = nca.shape, IndividualId = parse.(Int, string.(nca.id)),
                                Cmax = nca.cmax, AUC_tEnd = nca.auc)
end
CSV.write(joinpath(OUT, "step9_parameters.csv"), params)
CSV.write(joinpath(OUT, "step9_predicted_endpoints.csv"),
          reduce(vcat, [insertcols(d, 1, :method => m) for (m, d) in predicted]))
sources = ["PBPK" => pbpk; [m => predicted[m] for m in params.method]; "M5_all27" => m5]

# ---- 2. GMR of every formulation: predicted vs PBPK ----------------------------------
gmr = DataFrame()
for f in eachrow(FORMS), ep in ENDPOINTS
    vals = [Symbol("GMR_", m) => be_test(subject_values(pk, f.T50, f.shape, ep),
                                         subject_values(pk, REF.T50, REF.shape, ep)).GMR for (m, pk) in sources]
    push!(gmr, merge((T50 = f.T50, shape = f.shape, endpoint = ep), (; vals...)))
end
CSV.write(joinpath(OUT, "step9_gmr.csv"), gmr)

# ---- 3. BE decision per formulation (power >= 80%), as in step 5 ----------------------
dec = DataFrame()
for f in eachrow(FORMS), (m, pk) in sources
    test(ep) = subject_values(pk, f.T50, f.shape, ep)
    ref(ep)  = subject_values(pk, REF.T50, REF.shape, ep)
    push!(dec, (source = m, T50 = f.T50, shape = f.shape,
                ABE_joint = pass_rate_joint(test("Cmax"), ref("Cmax"), test("AUC_tEnd"), ref("AUC_tEnd");
                                            n_subj = N_ABE, n_trial = N_TRIAL_POWER,
                                            cvC = WSV_CV["Cmax"], cvA = WSV_CV["AUC_tEnd"]) >= POWER_TARGET,
                RSABE_Cmax = pass_rate_scaled(test("Cmax"), ref("Cmax"); n_subj = N_SCALED,
                                              n_trial = N_TRIAL_POWER, method = :fda_rsabe) >= POWER_TARGET,
                ABEL_Cmax  = pass_rate_scaled(test("Cmax"), ref("Cmax"); n_subj = N_SCALED,
                                              n_trial = N_TRIAL_POWER, method = :ema_abel) >= POWER_TARGET))
end
CSV.write(joinpath(OUT, "step9_decisions.csv"), dec)

# ---- 4. boundaries on each axis and the TRUE (PBPK) type I error there -----------------
bnd = DataFrame()
for ax in AXES, (m, pk) in sources, side in (:low, :high)
    G = grid_of(pk, ax)
    hits = [first_crossing(G, [gmr_at(pk, ax, ep, x) for x in G], ax.ref, side) for ep in ENDPOINTS]
    ok = [h for h in hits if !ismissing(h.x)]
    b = isempty(ok) ? missing : ok[argmin([abs(h.x - ax.ref) for h in ok])].x      # joint = first reached
    push!(bnd, (axis = ax.name, source = m, side = String(side), boundary = b,
                true_type1_n54 = ismissing(b) ? missing : pass_rate_at(pbpk, ax, "Joint", b, 54, N_TRIAL_TYPE1));
          promote = true)
end
CSV.write(joinpath(OUT, "step9_boundaries.csv"), bnd)

# ---- figure: Cmax GMR along each axis ---------------------------------------------------
COLOR = Dict("PBPK" => :black, "M5_all27" => :gray60,
             Pair.(params.method, [:royalblue, :dodgerblue, :navy, :seagreen, :limegreen, :darkgreen])...)
fig = Figure(size = (1100, 850))
for (i, ep) in enumerate(ENDPOINTS), (j, ax) in enumerate(AXES)
    a = Axis(fig[i, j]; title = "$ep GMR along the $(ax.name) axis", xlabel = ax.name,
             ylabel = "GMR vs reference", xscale = ax.log_x ? log10 : identity)
    hlines!(a, [BE_LOWER, BE_UPPER]; color = :gray70, linestyle = :dash)
    d = sort(gmr[(gmr.endpoint .== ep) .& (gmr[!, ax.fixed] .== ax.fixed_value), :], ax.col)
    for (m, _) in sources
        scatterlines!(a, d[!, ax.col], d[!, Symbol("GMR_", m)]; color = COLOR[m], label = m)
    end
    i == 1 && j == 1 && axislegend(a; position = :rt)
end
save(joinpath(OUT, "fig_step9_gmr.png"), fig)

println("Trained models:"); show(params; allcols = true); println()
println("\nGMR, predicted vs PBPK:"); show(gmr; allrows = true, allcols = true); println()
println("\nBE decisions (power >= 80%) vs PBPK:")
for (m, _) in sources[2:end], col in (:ABE_joint, :RSABE_Cmax, :ABEL_Cmax)
    a = dec[dec.source .== "PBPK", col]; b = dec[dec.source .== m, col]
    println(rpad(m, 14), rpad(col, 11), " BE: ", sum(b), " (PBPK ", sum(a), "), same: ", sum(a .== b),
            "/27, too lenient: ", sum(b .& .!a), ", too strict: ", sum(.!b .& a))
end
println("\nBoundaries and true (PBPK) type I error there:"); show(bnd; allrows = true, allcols = true); println()
