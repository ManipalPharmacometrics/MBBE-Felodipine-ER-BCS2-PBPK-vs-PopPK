# =============================================================================
# poppk_models_felo_structure.jl — structural variants of M3 aimed at the early
# Tmax (2026-10-24). All in the @delay form with the release lag,
# as src/poppk_models_felo_delay.jl (truncation TAU0 = 1e-4 h, see there).
#
# M3 (Weibull -> ka -> 2-cmt) rises correctly but peaks 0.4–0.7 h early and
# falls faster than PK-Sim at some arms (worst at shape 1.5): one release
# timescale + one ka cannot give a fast phase and a prolonged late tail.
#
#   poppk_m5  dual Weibull (FF*dd1 + (1-FF)*dd2) -> ka -> 2-cmt
#   poppk_m6  Weibull -> 3 transit compartments (ktr) -> 2-cmt (replaces ka)
#   poppk_m3d with tvlag estimated = M3L (no new model; fitted in 06)
#
# IIV as M3: one release eta (on both T50 in M5), CL, Vc. Units h, nmol, L, nmol/L.
# =============================================================================

poppk_m5 = @model begin
    @param begin
        tvt50a  ∈ RealDomain(lower = 0.05, init = 2.0)     # fast release component
        tvga    ∈ RealDomain(lower = 0.1,  init = 0.8)
        tvt50b  ∈ RealDomain(lower = 0.05, init = 10.0)    # slow release component
        tvgb    ∈ RealDomain(lower = 0.1,  init = 0.8)
        tvff    ∈ RealDomain(lower = 0.0, upper = 1.0, init = 0.5)   # fraction in the fast component
        tvka    ∈ RealDomain(lower = 0.01, init = 1.0)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 200.0)
        tvq     ∈ RealDomain(lower = 0.1,  init = 150.0)
        tvvp    ∈ RealDomain(lower = 1.0,  init = 1000.0)
        Ω       ∈ PDiagDomain(init = [0.05, 0.09, 0.09])
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        T50a = tvt50a * exp(η[1]); T50b = tvt50b * exp(η[1])
        CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Ga = tvga; Gb = tvgb; FF = tvff; Ka = tvka; Q = tvq; Vp = tvvp
        Aa = T50a / log(2)^(1 / Ga); Ab = T50b / log(2)^(1 / Gb)
        Tlag = tvlag; Tau0 = 1e-4                # truncation, see poppk_models_felo_delay.jl
        dd1  = @delay(Tlag + truncated(Weibull(Ga, Aa), Tau0, Inf), Gut)
        dd2  = @delay(Tlag + truncated(Weibull(Gb, Ab), Tau0, Inf), Gut)
    end
    @dosecontrol begin
        bioav = (; Gut = 0.0,)
    end
    @dynamics begin
        Gut'     = FF * dd1 + (1 - FF) * dd2 - Ka * Gut
        Central' = Ka * Gut - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  = Q / Vc * Central - Q / Vp * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end

poppk_m6 = @model begin
    @param begin
        tvt50   ∈ RealDomain(lower = 0.05, init = 4.0)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.8)
        tvktr   ∈ RealDomain(lower = 0.01, init = 3.0)     # transit rate (3 compartments)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 200.0)
        tvq     ∈ RealDomain(lower = 0.1,  init = 150.0)
        tvvp    ∈ RealDomain(lower = 1.0,  init = 1000.0)
        Ω       ∈ PDiagDomain(init = [0.05, 0.09, 0.09])
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        T50 = tvt50 * exp(η[1]); CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Gamma = tvgamma; Ktr = tvktr; Q = tvq; Vp = tvvp
        Alpha = T50 / log(2)^(1 / Gamma)
        Tlag  = tvlag; Tau0 = 1e-4
        dd    = @delay(Tlag + truncated(Weibull(Gamma, Alpha), Tau0, Inf), Tr1)
    end
    @dosecontrol begin
        bioav = (; Tr1 = 0.0,)
    end
    @dynamics begin
        Tr1'     = dd - Ktr * Tr1
        Tr2'     = Ktr * Tr1 - Ktr * Tr2
        Tr3'     = Ktr * Tr2 - Ktr * Tr3
        Central' = Ktr * Tr3 - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  = Q / Vc * Central - Q / Vp * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
