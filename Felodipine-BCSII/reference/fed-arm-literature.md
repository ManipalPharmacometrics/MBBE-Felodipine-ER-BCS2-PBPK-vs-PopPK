# Felodipine ER fed arm (gate #3): literature search, 2026-10-16

The goal is fed Plendil ER 10 mg data to validate a PK-Sim fed reference (scope §3.6, criterion 3: "food-effect direction and rough magnitude reproduced").
MW 384.26 is used for the unit conversions (1 ng/mL = 2.602 nM).

## Summary-level data (Plendil arm, single dose 10 mg, high-fat or high-carbohydrate meal)

The pdfs are in `literature-fed/`. The figures below are the **reference (Plendil) arm**, geometric means.

| Source | n | Fasted Cmax | Fed Cmax | Cmax fed/fasted | Fasted AUC∞ | Fed AUC∞ | AUC fed/fasted | Tmax fasted → fed |
|---|---|---|---|---|---|---|---|---|
| Sandoz CA monograph (3 BE studies: fasting, fed, SS; healthy males) | not stated | 2.45 ng/mL (6.38 nM) | 6.80 ng/mL (17.7 nM) | **2.77** | 37.2 ng·h/mL | 40.7 | **1.09** | median 4.5 → 4.5 h |
| Apotex CA monograph (separate studies, different subjects) | 49 fasted / 26 fed | 4.27 ng/mL (11.1 nM) | 15.1 ng/mL (39.4 nM) | 3.54* | 69.6 | 102.4 | 1.47* | mean 5.1 → 4.9 h |
| Plendil US label | — | — | — | "≈ +60%" | — | — | "unchanged" | — |
| Plendil / generic CA labels | — | — | — | "1.5–2 fold" | — | — | "not influenced" | — |
| Weitschies 2005 (J Control Release 108:375; numbers via Kim 2014 Table 3) | 6, crossover | 16.7 nM (arithmetic mean) | 27.3 nM | 1.63 | AUCall 53.9 nM·h | 54.3 | ≈1.0 (window < 8 h) | 4.1 → 4.5 h; tlag 0.4 → 1.6 h |

\* The Apotex fasted and fed studies used different subjects, and that population has about 2× the exposure in both arms. The cross-study ratio is confounded, so it is not used as a target.

**Model check (fasted, Pop_ref, GM of 100 subjects):** Cmax 5.16 nM = 1.98 ng/mL, ratio **0.81** vs Sandoz Plendil (consistent with the accepted 0.79 bias). AUC∞ 86.1 nM·h = 33.1 ng·h/mL, ratio **0.89**. tmax 2.7 h vs an observed median of 4.5 h: the model peaks early.

**Proposed fed validation target:** Cmax fed/fasted ratio between **1.6 and 2.8**, AUC ratio **≈ 1.0–1.1**, tmax not earlier than fasted.
- Lower end: US label and Weitschies.
- Upper end: Sandoz, the only Plendil within-program comparison.
- The profile shape should show the delayed, high late peaks that Weitschies saw in individual subjects.

## Profile-level data (for digitising)

- **Weitschies 2005 / Kim 2014 (PLoS One 9:e108392, open access, PMC4182452):**
  - Figs 3 and 4 give individual fasted and fed profiles for 6 subjects.
  - The tablet is a magnetically labelled felodipine ER hydrophilic-matrix tablet (Abrahamsson is a co-author; probably the Plendil formulation, not confirmed).
  - The observation window is short (< 8 h, as Kim 2014 reports), so the profiles cover absorption only.
- **Abrahamsson 1993 (Pharm Res 10:709):**
  - Plendil ER 10 mg in 8 healthy males, fasting vs non-fasting, gamma scintigraphy plus plasma.
  - Gastric emptying was 0.6 h fasted vs 3.2 h fed; there was colonic absorption; erosion was not affected by food.
  - **Paywalled (Springer).** This is the best match to the actual product.
- **Bergstrand 2009 (CPT 86:77):** a model of the same MMM study (Markov GI transit). Food lowered gastric transition and raised fraction absorbed across the gut wall. Mechanistic support only.

## Dead ends

- Torrent fed BE NCT01630902 and Ranbaxy fed BE vs Plendil NCT02311530 (n=80, high-fat breakfast 30 min pre-dose): **no results posted**.
- "Pharma Science Monitor" fed BE (India, Mylan reference): ResearchGate only, 403.
- The Plendil NDA 19-834 biopharm review is not online.
- Chinese felodipine SR fasting and fed studies: registrations only.
