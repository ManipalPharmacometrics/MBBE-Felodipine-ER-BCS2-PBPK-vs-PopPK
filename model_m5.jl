# =============================================================================
# model_m5.jl — the final felodipine PopPK model (M5), fitted per formulation.
#
#   Release: two Weibull components (fast "a", slow "b"); a fraction FF of the
#   dose follows the fast one. Both start after a fixed lag (0.3113 h).
#       fraction released (excluding the lag) =
#           FF * (1 - 2^-(t/T50a)^Ga) + (1 - FF) * (1 - 2^-(t/T50b)^Gb)
#   Written with Pumas @delay: each Weibull is the delay distribution of the dose
#   into the Gut compartment (Pumas `Weibull(shape, scale)`, scale = T50/ln2^(1/G)).
#   The dose itself gets bioavailability 0 into Gut; it only enters via @delay.
#
#       Gut'     = FF * input_a + (1 - FF) * input_b - Ka * Gut
#       Central' = Ka * Gut - (CL + Q)/Vc * Central + Q/Vp * Periph
#       Periph'  = Q/Vc * Central - Q/Vp * Periph
#       Cp       = Central / Vc
#
#   Between-subject variability: one eta on release (both T50s), CL, Vc.
#   Residual error: proportional.  Units: h, nmol, L, nmol/L (nM).
#   truncated(Weibull, 1e-4, Inf): avoids an infinite density at t = 0 when the
#   shape is < 1.
# =============================================================================

poppk_m5 = @model begin
    @param begin
        tvt50a ∈ RealDomain(lower = 0.05, init = 2.0)     # fast component T50 (h)
        tvga   ∈ RealDomain(lower = 0.1,  init = 0.8)     # fast component shape
        tvt50b ∈ RealDomain(lower = 0.05, init = 10.0)    # slow component T50 (h)
        tvgb   ∈ RealDomain(lower = 0.1,  init = 0.8)     # slow component shape
        tvff   ∈ RealDomain(lower = 0.0, upper = 1.0, init = 0.5)   # fraction in fast component
        tvka   ∈ RealDomain(lower = 0.01, init = 1.0)     # absorption rate (1/h)
        tvcl   ∈ RealDomain(lower = 1.0,  init = 300.0)   # clearance (L/h)
        tvvc   ∈ RealDomain(lower = 1.0,  init = 200.0)   # central volume (L)
        tvq    ∈ RealDomain(lower = 0.1,  init = 150.0)   # intercompartmental CL (L/h)
        tvvp   ∈ RealDomain(lower = 1.0,  init = 1000.0)  # peripheral volume (L)
        Ω      ∈ PDiagDomain(init = [0.05, 0.09, 0.09])   # IIV: release, CL, Vc
        tvlag  ∈ RealDomain(lower = 0.0,  init = 0.3113)  # release lag (h), held fixed
        σ_p    ∈ RealDomain(lower = 0.0,  init = 0.15)    # proportional error SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @pre begin
        T50a = tvt50a * exp(η[1]); T50b = tvt50b * exp(η[1])
        CL = tvcl * exp(η[2]); Vc = tvvc * exp(η[3])
        Ga = tvga; Gb = tvgb; FF = tvff; Ka = tvka; Q = tvq; Vp = tvvp
        Aa = T50a / log(2)^(1 / Ga)                # Weibull scale from T50
        Ab = T50b / log(2)^(1 / Gb)
        Tlag = tvlag; Tau0 = 1e-4
        dd1 = @delay(Tlag + truncated(Weibull(Ga, Aa), Tau0, Inf), Gut)
        dd2 = @delay(Tlag + truncated(Weibull(Gb, Ab), Tau0, Inf), Gut)
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
