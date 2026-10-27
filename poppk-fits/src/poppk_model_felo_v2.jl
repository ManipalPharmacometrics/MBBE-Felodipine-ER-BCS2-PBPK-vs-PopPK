# =============================================================================
# poppk_model_felo_v2.jl — candidate improved felodipine PopPK models (2026-10-23).
#
# Why: the 1-cmt Weibull baseline (poppk_model_felo_1cmt.jl) cannot produce the
# PBPK profile's fast rise + long plateau; it falls into a disposition-limited
# basin (release T50 ~0.6 h, Vc ~4000 L) that over-predicts Cmax ~10% and
# switches basin between arms (notes §9.27).
#
# Structure: Weibull release (parameterised in T50 so it can be compared with
# the PK-Sim T50) -> gut depot with first-order ka -> 2-cmt disposition.
#
#   poppk_felo_wka2   one formulation arm (per-arm fit, as the baseline)
#   poppk_felo_pair   reference + test in the same subjects (crossover-like):
#                     shared disposition + subject etas, formulation-specific
#                     T50 and shape. Period 2 is dosed at TP2 into Depot2.
#
# Units: h, nmol, L, nmol/L. Release lag held at the PK-Sim value 0.3113 h.
# The hazard clock is floored at 0.01 h: for shape < 1 the Weibull hazard is
# singular at the lag and an unfloored start gives an Inf likelihood.
# =============================================================================

const TP2 = 200.0   # period-2 dose time in the paired dataset (washout >> 5 t½)


poppk_felo_wka2 = @model begin
    @param begin
        tvt50   ∈ RealDomain(lower = 0.05, init = 4.0)
        tvgamma ∈ RealDomain(lower = 0.1,  init = 0.8)
        tvka    ∈ RealDomain(lower = 0.01, init = 1.0)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 200.0)
        tvq     ∈ RealDomain(lower = 0.1,  init = 150.0)
        tvvp    ∈ RealDomain(lower = 1.0,  init = 1000.0)
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        Ω       ∈ PDiagDomain(init = [0.05, 0.09, 0.09])
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        T50 = tvt50 * exp(η[1]); CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Gamma = tvgamma; Ka = tvka; Q = tvq; Vp = tvvp; Tlag = tvlag
        Alpha = T50 / log(2)^(1 / Gamma)          # Weibull scale from T50
    end
    @dynamics begin
        Depot'   = -ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 0.01)/Alpha)^(Gamma - 1), 0.0) * Depot
        Gut'     =  ifelse(t > Tlag, (Gamma/Alpha) * (max(t - Tlag, 0.01)/Alpha)^(Gamma - 1), 0.0) * Depot - Ka * Gut
        Central' =  Ka * Gut - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  =  Q / Vc * Central - Q / Vp * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end

poppk_felo_pair = @model begin
    @param begin
        t50R    ∈ RealDomain(lower = 0.05, init = 4.0)
        gR      ∈ RealDomain(lower = 0.1,  init = 0.8)
        t50T    ∈ RealDomain(lower = 0.05, init = 4.0)
        gT      ∈ RealDomain(lower = 0.1,  init = 0.8)
        tvka    ∈ RealDomain(lower = 0.01, init = 1.0)
        tvcl    ∈ RealDomain(lower = 1.0,  init = 300.0)
        tvvc    ∈ RealDomain(lower = 1.0,  init = 200.0)
        tvq     ∈ RealDomain(lower = 0.1,  init = 150.0)
        tvvp    ∈ RealDomain(lower = 1.0,  init = 1000.0)
        tvlag   ∈ RealDomain(lower = 0.0,  init = 0.3113)
        Ω       ∈ PDiagDomain(init = [0.05, 0.09, 0.09])
        σ_p     ∈ RealDomain(lower = 0.0,  init = 0.15)
    end
    @random begin
        η ~ MvNormal(Ω)
    end
    @pre begin
        T50R = t50R * exp(η[1]); T50T = t50T * exp(η[1])   # subject release eta shared (GI transit)
        CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Ka = tvka; Q = tvq; Vp = tvvp; Tlag = tvlag; GR = gR; GT = gT
        AR = T50R / log(2)^(1 / GR); AT = T50T / log(2)^(1 / GT); Tp2 = TP2
    end
    @dynamics begin
        Depot1'  = -ifelse(t > Tlag, (GR/AR) * (max(t - Tlag, 0.01)/AR)^(GR - 1), 0.0) * Depot1
        Depot2'  = -ifelse(t > Tp2 + Tlag, (GT/AT) * (max(t - Tp2 - Tlag, 0.01)/AT)^(GT - 1), 0.0) * Depot2
        Gut'     =  ifelse(t > Tlag, (GR/AR) * (max(t - Tlag, 0.01)/AR)^(GR - 1), 0.0) * Depot1 +
                    ifelse(t > Tp2 + Tlag, (GT/AT) * (max(t - Tp2 - Tlag, 0.01)/AT)^(GT - 1), 0.0) * Depot2 - Ka * Gut
        Central' =  Ka * Gut - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  =  Q / Vc * Central - Q / Vp * Periph
    end
    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end

# Starting points: fast-absorption / slow-release (release-limited) and the
# reverse (absorption-limited); fits keep the higher log-likelihood.
const INITS_V2 = [
    "release_limited"    => (tvka = 2.0, tvcl = 300.0, tvvc = 200.0, tvq = 200.0, tvvp = 1500.0),
    "absorption_limited" => (tvka = 0.3, tvcl = 300.0, tvvc = 500.0, tvq = 100.0, tvvp = 1000.0)]
