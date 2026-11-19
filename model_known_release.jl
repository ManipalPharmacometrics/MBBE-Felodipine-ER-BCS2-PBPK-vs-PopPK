# =============================================================================
# model_known_release.jl — "known release + regional gut availability" PopPK
# model for felodipine (step 9).
#
# The release is NOT estimated: it is each formulation's own PK-Sim dissolution
# curve, given as data (covariates T50 and shape; lag 0.3113 h as in PK-Sim):
#     fraction released = 1 - 2^-(s/T50)^shape,   s = time since the lag
#     release rate h(s) = ln2 * shape * s^(shape-1) / T50^shape
# What the model estimates is what happens after release, the same for every
# formulation:
#   - first-order absorption (ka);
#   - gut-wall availability that changes as the tablet moves down the gut.
#     Felodipine is a CYP3A4 substrate; gut-wall CYP3A4 is high in the small
#     intestine and low in the colon, so drug released late escapes more of it:
#         r(t) = 1 + (R - 1) * (1 - exp(-t / tau))
#     r = 1 early (small intestine), r -> R late (colon). R > 1 means drug
#     absorbed in the colon is more available; tau = time scale of the move.
#     tau is FIXED in step 9: estimated freely it ran to ~0.2 h, where R just
#     scales availability overall and trades off against CL and the volumes
#     (oral data cannot separate them);
#   - 2-compartment disposition (CL, Vc, Q, Vp). CL and volumes are apparent
#     (relative to small-intestine availability), as with any oral-only data.
#
#   Tablet'  = -h(t - lag) * Tablet                          dose enters at the lag
#   Gut'     =  h(t - lag) * Tablet - ka * Gut
#   Central' =  r(t) * ka * Gut - (CL + Q)/Vc * Central + Q/Vp * Periph
#   Periph'  =  Q/Vc * Central - Q/Vp * Periph
#
# IIV: tau, CL, Vc. Proportional error. Units: h, nmol, L, nM.
# To predict an untested formulation: set its T50 and shape, keep everything else.
# =============================================================================

poppk_known_release = @model begin
    @param begin
        tvka  ∈ RealDomain(lower = 0.01, init = 3.0)     # absorption rate (1/h)
        tvR   ∈ RealDomain(lower = 0.0,  init = 1.5)     # colon / small-intestine availability (0 = none)
        tvtau ∈ RealDomain(lower = 0.1,  init = 5.0)     # time scale of the move to the colon (h)
        tvcl  ∈ RealDomain(lower = 1.0,  init = 300.0)   # clearance (L/h)
        tvvc  ∈ RealDomain(lower = 1.0,  init = 70.0)    # central volume (L)
        tvq   ∈ RealDomain(lower = 0.1,  init = 20.0)    # intercompartmental CL (L/h)
        tvvp  ∈ RealDomain(lower = 1.0,  init = 150.0)   # peripheral volume (L)
        Ω     ∈ PDiagDomain(init = [0.09, 0.09, 0.09])   # IIV: tau, CL, Vc
        σ_p   ∈ RealDomain(lower = 0.0,  init = 0.15)    # proportional error SD
    end

    @random begin
        η ~ MvNormal(Ω)
    end

    @covariates T50 shape

    @pre begin
        Ka  = tvka
        R   = tvR
        Tau = tvtau * exp(η[1])
        CL  = tvcl * exp(η[2])
        Vc  = tvvc * exp(η[3])
        Q   = tvq
        Vp  = tvvp
        B   = shape
        Rel = log(2) * shape / T50^shape          # known release: h(s) = Rel * s^(shape-1)
    end

    @dosecontrol begin
        lags = (; Tablet = 0.3113,)                  # PK-Sim tablet lag (18.68 min)
    end

    @dynamics begin
        Tablet'  = -Rel * max(t - 0.3113, 1e-8)^(B - 1) * Tablet
        Gut'     =  Rel * max(t - 0.3113, 1e-8)^(B - 1) * Tablet - Ka * Gut
        Central' =  (1 + (R - 1) * (1 - exp(-t / Tau))) * Ka * Gut - (CL + Q) / Vc * Central + Q / Vp * Periph
        Periph'  =  Q / Vc * Central - Q / Vp * Periph
    end

    @derived begin
        cp := @. Central / Vc
        dv ~ @. Normal(cp, abs(cp) * σ_p + 1e-6)
    end
end
