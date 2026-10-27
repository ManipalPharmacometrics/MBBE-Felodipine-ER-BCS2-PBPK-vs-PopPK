# =============================================================================
# poppk_models_felo_delay.jl — M1/M2/M3 rewritten in the Pumas @delay form
# (dose into the receiving
# compartment with bioav 0, input = @delay(Weibull(shape, scale), cmt))
# (2026-10-24).
#
# Purpose: check whether the fitted model results stay the same when the
# hand-coded Weibull hazard ODE (M1–M3) is replaced by Pumas @delay.
#
# Release lag: input = @delay(Tlag + truncated(Weibull(Gamma, Alpha), TAU0, Inf)).
#   - @dosecontrol lags are ignored by @delay (prototype: identical output).
#   - An untruncated Tlag + Weibull, or a dose record at t = Tlag, returns Inf
#     for shape < 1 (density singular at the dose time).
#   - Truncating at TAU0 = 1e-4 h is finite for every shape and removes ≤0.4% of
#     the release mass (shape 0.5). It plays the role of the ODE models' clock
#     floor (1e-8 h in M1/M2, 0.01 h in M3), and is closer to the exact Weibull
#     than the 0.01 h floor (prototype: within ~1% of the ODE from 2 h on).
#   tvlag is held constant in the fits: 0.3113 h (PK-Sim value, as M1–M3) or 0
#   for the no-lag variant (the plain @delay form has no lag).
# Prototypes: session scratchpad proto_delay.jl, proto_nolag.jl, proto_lag.jl.
#
# Units: h, nmol, L, nmol/L. IIV and residual error as in M1–M3.
# =============================================================================

# M1d: 1-cmt, Weibull input in scale form (as M1)
poppk_m1d = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.05, init = 7.0)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.9)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 800.0)
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        Alpha = tvalpha * exp(η[1]); CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Gamma = tvgamma
        Tlag  = tvlag; Tau0 = 1e-4                # truncation, see header
        dd    = @delay(Tlag + truncated(Weibull(Gamma, Alpha), Tau0, Inf), Central)
    end
    @dosecontrol begin
        bioav = (; Central = 0.0,)
    end
    @dynamics begin
        Central' = dd - (CL / Vc) * Central
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end

# M2d: 2-cmt, Weibull input in scale form (as M2)
poppk_m2d = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.05, init = 7.0)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.9)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 250.0)
        tvq     ∈ RealDomain(lower = 0.1,  init = 200.0)
        tvvp    ∈ RealDomain(lower = 1.0,  init = 1500.0)
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        Alpha = tvalpha * exp(η[1]); CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Gamma = tvgamma; Q = tvq; Vp = tvvp
        Tlag  = tvlag; Tau0 = 1e-4                # truncation, see header
        dd    = @delay(Tlag + truncated(Weibull(Gamma, Alpha), Tau0, Inf), Central)
    end
    @dosecontrol begin
        bioav = (; Central = 0.0,)
    end
    @dynamics begin
        Central' = dd - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  = Q / Vc * Central - Q / Vp * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end

# M3d: Weibull input (T50 form, as M3) into the gut depot -> ka -> 2-cmt
poppk_m3d = @model begin
    @param begin
        tvt50   ∈ RealDomain(lower = 0.05, init = 4.0)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.8)
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
        T50 = tvt50 * exp(η[1]); CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Gamma = tvgamma; Ka = tvka; Q = tvq; Vp = tvvp
        Alpha = T50 / log(2)^(1 / Gamma)
        Tlag  = tvlag; Tau0 = 1e-4                # truncation, see header
        dd    = @delay(Tlag + truncated(Weibull(Gamma, Alpha), Tau0, Inf), Gut)
    end
    @dosecontrol begin
        bioav = (; Gut = 0.0,)
    end
    @dynamics begin
        Gut'     = dd - Ka * Gut
        Central' = Ka * Gut - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  = Q / Vc * Central - Q / Vp * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
