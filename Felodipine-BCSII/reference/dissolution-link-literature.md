# Felodipine ER dissolution link (gate #4): literature search, 2026-10-16

Goal: make the reference in-vivo Weibull (T50 4.78 h, shape 0.764, lag 18.68 min) traceable to real Plendil ER in-vitro dissolution (scope criterion 4).

## USP monograph (USP35, Felodipine ER Tablets)

| Test | Conditions | 2 h | 6 h | 10 h |
|---|---|---|---|---|
| Test 1 | App. 2, 50 rpm, 500 mL pH 6.5 phosphate + 1% SLS, stationary basket | 10–30% | 42–68% | ≥ 75% |
| Test 3 | App. 2, 50 rpm, 500 mL pH 6.5 phosphate + 1% SLS, stationary basket | 10–30% | 50–80% | ≥ 80% |

Test 2 (App. 1, 100 rpm, 1% polysorbate 80) checks at 1, 4 and 8 h instead: 5–30%, 45–70%, ≥ 80%.

**Reference in-vivo Weibull, read as if it were an in-vitro curve:** 26.9%, 54.7%, **69.6%** at 2, 6 and 10 h.
- It passes at 2 h and 6 h but **fails the 10 h limit**.
- Real Plendil passes USP, so in-vitro release is faster and closer to zero-order late on than the in-vivo release the plasma data imply.
- This fits Abrahamsson 1993: erosion is slower in the distal gut. It also fits Javed 2022: Plendil in-vitro release is zero-order.

**Conclusion:** in vitro and in vivo differ in **shape**, not only in timescale. A 1:1 mapping is wrong, so gate #4 needs a fitted in-vitro↔in-vivo relation (Levy plot).

If the curves are read 1:1, the (T50, shape) region that meets Test 1 is:

| Shape | T50 range |
|---|---|
| ≤ 0.764 | none |
| 0.8 | 3.9–4.1 h |
| 1.0 | 3.5–4.8 h |
| 1.2 | 3.8–5.4 h |
| 1.5 | 4.1–5.9 h |

This is only partly inside the BE band. It is not used, for the reason above.

## Sources

- **Abrahamsson 1994, Pharm Res 11:1093 (PMID 7971707). KEY, paywalled.**
  - Three slightly different felodipine ER 10 mg gel-matrix tablets (the Plendil platform), tested in vitro (SLS, Tween, CTAB media) and in a bioavailability study with an oral solution.
  - In-vitro vs in-vivo dissolution: r² 0.94–0.97 by moment analysis. SLS gives a much steeper in-vitro/in-vivo slope.
  - This gives both the IVIVC and a real three-formulation in-vivo check of the safe space.
- **Abrahamsson 1993, Pharm Res 10:709.** In-vivo erosion by scintigraphy, slower distally. Paywalled (also listed in the fed note).
- **Vetchy 2007, Medicina (Kaunas) 43:326 (MDPI, open access).** Plendil ER 5 and 10 mg vs Presid, 24 h profiles at pH 1.2, 4.5 and 6.5. The automated download was blocked (Access Denied), so it needs to be downloaded in a browser.
- **Javed 2022, Curr Pharm Des 28:2010.** Plendil 5 mg profiles, standard USP method vs 1% Tween 20; zero-order release. Paywalled.
- **Petersen 2003, Arzneimittelforschung 53:40.** Modip 5 mg (the German brand name for Plendil) vs 7 generics, 7 h, linear release. Paywalled.
- **Wingstrand 1990** (bioavailability from felodipine ER tablets with different dissolution properties, Int J Pharm). Not found in PubMed; probably Int J Pharm 60:151.

## Update 2026-10-16 (pm): papers read

PDFs (not included): Abrahamsson 1994 (a scan, read as images) and Vetchy 2007. **Abrahamsson 1993 was not uploaded.**
Digitised data: `digitised/` (see its README for accuracy).

### Abrahamsson 1994: what it gives

- **Tablets A, B and C:** identical gel-matrix (HPMC) felodipine ER 10 mg tablets that differ only in HPMC viscosity. They were tested in a 4-way crossover against an oral solution, n=16 healthy males, 24 h.
- **Table III (exact):**

  | Tablet | MDTvivo (h) | AUC (nmol·h/L) | Frel (%) |
  |---|---|---|---|
  | A | 3.3 | 88 | 88 |
  | B | 7.6 | 92 | 95 |
  | C | 11.5 | 86 | 87 |
  | Solution | — | 99 | 100 |

- **Fig 3, peaks of the mean curves:** A 9.5 nM at 3 h, B 5.0 nM at 5 h, C 2.3 nM at 5 h, solution 34 nM at 0.5 h. Ratios: A/B 1.9, B/C 2.2.
- **In-vitro/in-vivo correlation, MDTvivo = f(MDTvitro), Table IV.** The **1% SLS, 50 rpm** row is the USP Test 1 condition.

  | Medium | 100 rpm | 50 rpm |
  |---|---|---|
  | 1% SLS | 8.1x − 14.8 | 2.9x − 7.4 |
  | 2% Tween | 2.1x − 0.5 | 1.1x + 0.1 |
  | 0.4% CTAB | 2.4x − 0.9 | 1.6x − 1.2 |

- **In-vitro release** is zero-order up to 90%, so MDTvitro is the time to 50% by linear regression.
- **Weibull fits to the Fig 4 in-vivo dissolution,** with lag fixed at 18.68 min and Fmax free:

  | Tablet | Fmax | T50 | Shape | Note |
  |---|---|---|---|---|
  | A | 0.80 | 1.70 h | 1.18 | |
  | B | 0.89 | 6.05 h | 0.94 | |
  | C | 1.00 (bound) | 19.75 h | 1.10 | Poorly identified: only 60% by 24 h |

  **The OSP Plendil reference (4.78 h, shape 0.764) sits between A and B, close to B.**

### Vetchy 2007

- Plendil ER 10 mg in 0.4% CTAB, paddle 100 rpm, **no stationary basket mentioned**.
- Rebuilt cumulative release at pH 6.5: 5 / 29 / 56 / 84 / 96 / 100 % at 1 / 4 / 7 / 12 / 18 / 24 h. MDTvitro 6.4–6.9 h.
- In Abrahamsson's equivalent medium (with the basket), tablets A, B and C give 1.9, 3.4 and 5.3 h. That would put Plendil slower than tablet C, which contradicts Plendil's in-vivo behaviour (tmax ≈ 4.5 h, Cmax ≈ 6.4 nM).
- **So Vetchy's set-up cannot be chained to the correlation**; the hydrodynamics differ without the basket. It only confirms near-zero-order release over 24 h.

### Gate #4 verdict: qualitatively consistent, not a quantitative IVIVC

- **USP Test 1 window → MDTvivo via 2.9x − 7.4:**
  - Fast edge: MDTvitro 4.1 h → 4.5 h in vivo.
  - Midpoint: 5.4 h → 8.3 h.
  - Slow edge: 7.0 h → 12.9 h.
- **Reference Weibull MDT** (infinite-time moment) is 9.4 h, inside the window and near its middle.
  - Caveat: Weibull MDTs for A, B and C (2.4, 10.2, 26 h) differ from Abrahamsson's tabulated values (3.3, 7.6, 11.5 h), so the moment definitions don't match. The comparison is approximate.
- **The USP-implied window is much wider than the BE band** (MDT 8.6–10.2 h at the reference shape).
  - The dissolution spec cannot by itself guarantee BE; this is a useful conclusion.
  - MDT alone doesn't decide BE either: across the band, MDT runs from 7 to 17 h depending on shape.

### Stronger check proposed: external validation with tablets A, B and C

Run 3 PK-Sim clones of Pop_ref with the fitted Weibulls (A 102.0 min / shape 1.18; B 363.0 min / 0.94; C 1185 min / 1.10). Compare simulated Cmax ratios A:B:C (observed ≈ 1.9 and 2.2) and AUC (observed ≈ equal, Frel 0.87–0.95).
- This tests the model's Cmax sensitivity to release rate, which is the thing that sets the width of the safe space.
- It does not depend on matching in-vitro methods.
