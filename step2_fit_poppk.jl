# =============================================================================
# STEP 2 — fit the M5 PopPK model to each of the 27 formulations (needs Pumas).
#
# In : results/step1_dataset.csv, model_m5.jl
# Out: results/fits/<T50>_<shape>/parameters.csv        estimates and -2 log-likelihood
#      results/fits/<T50>_<shape>/standard_errors.csv   estimates with SE and 95% CI
#      results/fits/<T50>_<shape>/ipred.csv             individual predictions, 0–96 h
#                                                       every 0.05 h (used in step 3;
#                                                       large, not committed)
#      results/step2_parameters.csv                     all formulations, once all are done
#
# Each formulation = 100 subjects, fitted on its own with FOCE. Lag fixed at
# 0.3113 h (the PK-Sim tablet lag). Starting values depend on the formulation's
# PK-Sim T50 (see `start`). One fit takes ~20–30 min on one core.
#
# Usage:  julia step2_fit_poppk.jl                 all formulations
#         julia step2_fit_poppk.jl 4.78_0.764 ...  only these (to run in parallel)
# Formulations already fitted (parameters.csv present) are skipped.
# =============================================================================

using Pumas, CSV, DataFrames
include("settings.jl")
include("model_m5.jl")

const PRED_TIMES = collect(0.0:0.05:96.0)

"Starting values (absorption-limited start) scaled by the formulation's T50."
start(T50) = (tvt50a = 0.3T50, tvga = 1.0, tvt50b = 2T50, tvgb = 1.0, tvff = 0.3, tvka = 0.5,
              tvcl = 300.0, tvvc = 500.0, tvq = 100.0, tvvp = 1000.0, tvlag = 0.3113)

data  = CSV.read(joinpath(RESULTS, "step1_dataset.csv"), DataFrame)
data.cmt .= 1                                   # all rows on Gut: the dose enters via @delay
forms = sort(unique(data[:, [:T50, :shape]]), [:shape, :T50])
name(f) = "$(f.T50)_$(f.shape)"
isempty(ARGS) || (forms = forms[in.(name.(eachrow(forms)), Ref(ARGS)), :])

for f in eachrow(forms)
    dir = joinpath(RESULTS, "fits", name(f))
    isfile(joinpath(dir, "parameters.csv")) && (println("skip ", name(f)); continue)
    mkpath(dir)

    arm = data[(data.T50 .== f.T50) .& (data.shape .== f.shape), :]
    pop = read_pumas(arm; id = :SUBJID, time = :time, amt = :amt, evid = :evid, cmt = :cmt,
                     observations = [:dv], covariates = [:T50, :shape, :IndividualId])
    @info "Fitting T50 = $(f.T50) h, shape = $(f.shape)"
    fitted = fit(poppk_m5, pop, merge(init_params(poppk_m5), start(f.T50)), FOCE();
                 constantcoef = (:tvlag,), optim_options = (show_trace = false,))

    est = coef(fitted)
    CSV.write(joinpath(dir, "standard_errors.csv"), DataFrame(coeftable(infer(fitted))))
    ipred = DataFrame()
    for (subj, pred) in zip(pop, predict(fitted; obstimes = PRED_TIMES))
        append!(ipred, DataFrame(T50 = f.T50, shape = f.shape,
                                 IndividualId = Int(subj.covariates(0.0).IndividualId),
                                 time = pred.time, ipred = pred.ipred.dv))
    end
    CSV.write(joinpath(dir, "ipred.csv"), ipred)
    CSV.write(joinpath(dir, "parameters.csv"),           # written last = formulation done
              DataFrame(T50 = f.T50, shape = f.shape, minus2LL = -2 * loglikelihood(fitted),
                        tvt50a = est.tvt50a, tvga = est.tvga, tvt50b = est.tvt50b, tvgb = est.tvgb,
                        tvff = est.tvff, tvka = est.tvka, tvcl = est.tvcl, tvvc = est.tvvc,
                        tvq = est.tvq, tvvp = est.tvvp, tvlag = est.tvlag,
                        omega_release = est.Ω[1, 1], omega_cl = est.Ω[2, 2], omega_vc = est.Ω[3, 3],
                        sigma_p = est.σ_p))
    println("done ", name(f)); flush(stdout)
end

# combined table once every formulation is fitted
all_dirs = [joinpath(RESULTS, "fits", d) for d in readdir(joinpath(RESULTS, "fits"))]
if count(d -> isfile(joinpath(d, "parameters.csv")), all_dirs) == 27
    params = reduce(vcat, [CSV.read(joinpath(d, "parameters.csv"), DataFrame) for d in all_dirs])
    CSV.write(joinpath(RESULTS, "step2_parameters.csv"), sort(params, [:shape, :T50]))
    println("Wrote results/step2_parameters.csv")
end
