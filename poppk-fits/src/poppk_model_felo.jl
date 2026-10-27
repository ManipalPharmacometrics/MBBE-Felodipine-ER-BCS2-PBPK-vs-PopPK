# =============================================================================
# poppk_model_felo.jl  — PopPK model for the felodipine ER PBPK -> PopPK
# decision-equivalence comparison (MBBE drug two).
#
# Same idea as the metformin baseline (metformin_poppk_weibull_export_edited/
# poppk_pumas/src/poppk_model.jl): the absorption input IS a Weibull in-vivo
# release process, so Alpha/Gamma are directly comparable to the PK-Sim
# T50/shape. Felodipine differs in disposition: it declines in 3 phases with a
# long terminal half-life (~25 h), so the disposition is 2-compartment.
#
#   Depot(t)  = Dose * exp(-((t - Tlag)/Alpha)^Gamma)         (Weibull release)
#   h(t)      = (Gamma/Alpha) * ((t - Tlag)/Alpha)^(Gamma-1)   (release hazard)
#   Central'  = h*Depot - (CL/Vc)*Central - (Q/Vc)*Central + (Q/Vp)*Periph
#   Periph'   = (Q/Vc)*Central - (Q/Vp)*Periph
#
# Alpha is the Weibull SCALE; the PK-Sim T50 is Alpha*(ln 2)^(1/Gamma).
# CL, Vc, Q, Vp are apparent (/F): felodipine F ≈ 0.15-0.2.
# Units: time h, amount nmol (10 mg = 26 024 nmol), volume L, conc nmol/L.
# IIV (diagonal) on Alpha, CL, Vc; proportional residual error.
# Tlag is held constant at the PK-Sim product lag (18.68 min = 0.3113 h).
# =============================================================================

poppk_model_felo = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.05, init = 7.0)    # Weibull release scale (h)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.9)    # Weibull shape
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)  # CL/F (L/h)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 800.0)  # Vc/F (L)
        tvq     ∈ RealDomain(lower = 1.0,  init = 400.0)  # Q/F (L/h)
        tvvp    ∈ RealDomain(lower = 1.0,  init = 4000.0) # Vp/F (L)
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113) # release lag (h), held constant
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])  # IIV: Alpha, CL, Vc
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)   # proportional residual SD
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        Alpha = tvalpha * exp(η[1])
        CL    = tvcl    * exp(η[2])
        Vc    = tvvc    * exp(η[3])
        Q     = tvq
        Vp    = tvvp
        Gamma = tvgamma
        Tlag  = tvlag
    end
    @dynamics begin
        Depot'   = -ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot
        Central' =  ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot -
                    (CL / Vc) * Central - (Q / Vc) * Central + (Q / Vp) * Periph
        Periph'  =  (Q / Vc) * Central - (Q / Vp) * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
