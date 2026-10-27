# Tests for be_utils.jl (2026-10-15 multi-drug refactor).
# Usage, from programs:  julia --project=. julia/test_be_utils.jl "$PWD/julia" <felodipine PK-Analyses dir>  (absolute paths)
using CSV, DataFrames, Statistics, Test
include(joinpath(ARGS[1], "be_utils.jl")); using .BEUtils
const FELO = ARGS[2]

@testset "sample size 2x2 ABE (PowerTOST sampleN.TOST, GMR 0.95, 80%)" begin
    @test sample_size_abe_2x2(0.20) == 20
    @test sample_size_abe_2x2(0.25) == 28
    @test sample_size_abe_2x2(0.30) == 40
end
for cv in (0.2434, 0.2525, 0.3571, 0.3651)
    println("n(2x2 ABE, CV=$cv, GMR 0.95, 80%) = ", sample_size_abe_2x2(cv),
            "   (GMR 0.90: ", sample_size_abe_2x2(cv; gmr = 0.90), ")")
end

@testset "parse_formulation" begin
    @test parse_formulation("Metformin_T50_3.5h_Simpop-PK-Analyses.csv", :metformin) |> x -> (x.T50 == 3.5 && ismissing(x.shape))
    @test parse_formulation("Pop_T50_2.39h-PK-Analyses.csv", :felodipine) == (T50 = 2.39, shape = 0.764)
    @test parse_formulation("Pop_b0.5-PK-Analyses.csv", :felodipine) == (T50 = 4.78, shape = 0.5)
    @test parse_formulation("Pop_b1.0-PK-Analyses.csv", :felodipine) == (T50 = 4.78, shape = 1.0)
    @test parse_formulation("Pop_T50_3.35h_b0.9-PK-Analyses.csv", :felodipine) == (T50 = 3.35, shape = 0.9)
    @test_throws ErrorException parse_formulation("Felodipine_T50_4.78h_b0764_x-PK-Analyses.csv", :felodipine)
end

@testset "wsv lookup" begin
    @test wsv_cv(:metformin) === WSV_CV
    @test wsv_cv(:metformin, :ir_sensitivity) === WSV_CV_IR
    @test wsv_cv(:felodipine, :fasted, "Cmax") == 0.3571
    @test wsv_cv(:felodipine, :fed, "AUC_inf") == 0.2216
end

pk  = load_all_formulations(FELO; drug = :felodipine)
println(combine(groupby(pk, [:T50, :shape]), nrow => :n,
        :Cmax => (x -> exp(mean(log.(x)))) => :gm_Cmax))
ref = reference_rows(pk, :felodipine)
tst = filter(r -> r.T50 == 2.39, pk)
be  = be_all_metrics(tst, ref)
println(be)
@testset "felodipine load" begin
    @test nrow(ref) == 100
    @test isapprox(be[be.metric .== "Cmax", :GMR][1], 1.488; atol = 0.001)   # python check
    @test isapprox(be[be.metric .== "AUC_inf", :GMR][1], 1.039; atol = 0.001)
end

cvC = wsv_cv(:felodipine, :fasted, "Cmax")
for nss in (sample_size_abe_2x2(cvC), 86)
    r1 = vbe_replicate(tst.Cmax, ref.Cmax; n_subj = nss, n_trial = 500, wsv_cv = cvC)
    println("ABE Cmax T50 2.39h n=$nss: pass ", r1.pass_rate, "  GMRmed ", round(r1.GMR_med, digits = 3))
end

@testset "scaled BE sanity" begin
    # ABEL cap: CVwR 50% -> 69.84-143.19%
    @test isapprox(exp(0.760 * sqrt(log(1.25))), 1.4319; atol = 1e-4)
    # self-comparison (T == R), CVwR 0.34, n 24: should pass most trials, mostly scaled
    s = vbe_replicate_scaled(ref.Cmax, ref.Cmax; n_subj = 24, n_trial = 2000, cv_wr = 0.34)
    a = vbe_replicate_scaled(ref.Cmax, ref.Cmax; n_subj = 24, n_trial = 2000, cv_wr = 0.34, method = :ema_abel)
    println("self RSABE n24: ", s); println("self ABEL  n24: ", a)
    @test s.pass_rate > 0.8 && a.pass_rate > 0.8
    @test 0.6 < s.scaled_rate < 1.0
    # low variability -> never scaled, falls back to ABE
    l = vbe_replicate_scaled(ref.Cmax, ref.Cmax; n_subj = 24, n_trial = 500, cv_wr = 0.10)
    @test l.scaled_rate == 0 && l.pass_rate > 0.99
    # true GMR 1.25 at CVwR 0.34: ABE ~ alpha, RSABE known to exceed it
    t125 = ref.Cmax .* 1.25
    abe = vbe_replicate(t125, ref.Cmax; n_subj = 48, n_trial = 2000, wsv_cv = 0.34)
    rs  = vbe_replicate_scaled(t125, ref.Cmax; n_subj = 24, n_trial = 2000, cv_wr = 0.34)
    println("GMR 1.25: ABE 2x2 n48 pass ", abe.pass_rate, " | RSABE n24 pass ", rs.pass_rate)
    @test abe.pass_rate < 0.07
end
