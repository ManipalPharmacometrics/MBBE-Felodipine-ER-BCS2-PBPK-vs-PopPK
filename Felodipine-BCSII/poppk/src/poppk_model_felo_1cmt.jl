# =============================================================================
# poppk_model_felo_1cmt.jl — BASELINE felodipine PopPK model (2026-10-16).
# 1-compartment disposition + Weibull in-vivo release input, structurally the
# same as the metformin baseline (metformin_poppk_weibull_export_edited/
# poppk_pumas/src/poppk_model.jl), with the felodipine release lag held at the
# PK-Sim product value (0.3113 h).
#
# Why not 2-compartment: on the reference arm both a disposition-like and an
# ER-like start drove Q or Vp to their bounds (Vp -> 1 L, or Q -> 5 L/h); the
# rich 0-96 h data do not identify a peripheral compartment
# (project notes §9.26). Two basins exist in 1-cmt as well: release-limited
# (flip-flop, slow Alpha) and disposition-limited (fast Alpha, large Vc). 02
# fits both starts per arm and keeps the higher log-likelihood.
#
# Units: h, nmol, L, nmol/L.  IIV on Alpha, CL, Vc.  Proportional error.
# =============================================================================

poppk_model_felo = @model begin
    @param begin
        tvalpha ∈ RealDomain(lower = 0.05, init = 7.0)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.9)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 800.0)
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        Ω       ∈ PDiagDomain(init = [0.09, 0.09, 0.09])
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        Alpha = tvalpha * exp(η[1])
        CL    = tvcl    * exp(η[2])
        Vc    = tvvc    * exp(η[3])
        Gamma = tvgamma
        Tlag  = tvlag
    end
    @dynamics begin
        Depot'   = -ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot
        Central' =  ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 1e-8)/Alpha)^(Gamma - 1), 0.0) * Depot - (CL / Vc) * Central
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end

# The two starting points (see header).
const INITS_FELO = [
    "release_limited"     => (tvalpha = 13.0, tvgamma = 1.0, tvcl = 300.0, tvvc = 250.0),
    "disposition_limited" => (tvalpha = 1.0,  tvgamma = 1.3, tvcl = 300.0, tvvc = 4000.0)]
