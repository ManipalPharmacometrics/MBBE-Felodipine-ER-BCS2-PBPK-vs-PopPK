// Felodipine ER (MBBE drug 2): full step-by-step work summary deck.
// One content spec -> two outputs: a .pptx (pptxgenjs) and a self-contained .html.
// Usage: NODE_PATH=<dir with pptxgenjs> node scripts/build_full_deck.js [outdir]
const fs = require("fs");
const path = require("path");
const pptxgen = require("pptxgenjs");

const ROOT = path.resolve(__dirname, "..");
const OUT = path.resolve(process.argv[2] || path.join(ROOT, "slides", "full-summary"));
const BASENAME = "Felodipine-ER-Full-Work-Summary";
fs.mkdirSync(OUT, { recursive: true });

// ---------- theme ----------
const THEME = {
  name: "Felodipine MBBE",
  headFontFace: "Cambria",
  bodyFontFace: "Calibri",
  colors: {
    dk1: "1B2430", lt1: "FFFFFF", dk2: "0F3D4C", lt2: "EDF3F4",
    accent1: "0F766E", accent2: "D97706", accent3: "64748B",
    accent4: "B42318", accent5: "15803D", accent6: "1D4ED8",
    hlink: "0F766E", folHlink: "64748B",
  },
};
const HEX = {
  text1: THEME.colors.dk1, text2: THEME.colors.dk2, background1: THEME.colors.lt1,
  background2: THEME.colors.lt2, accent1: THEME.colors.accent1, accent2: THEME.colors.accent2,
  accent3: THEME.colors.accent3, accent4: THEME.colors.accent4, accent5: THEME.colors.accent5,
  accent6: THEME.colors.accent6,
};
const W = 13.333, H = 7.5;
const CONTENT = { x: 0.6, y: 1.75, w: W - 1.2, h: 5.15 };
const GAP = 0.3;

// ---------- small helpers ----------
const img = (rel) => path.join(ROOT, rel);
function pngSize(file) {
  const b = fs.readFileSync(file);
  return { w: b.readUInt32BE(16), h: b.readUInt32BE(20) };
}
function esc(s) {
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}
// "**bold** rest" -> runs
function runs(s) {
  const out = [];
  String(s).split(/(\*\*[^*]+\*\*)/).forEach((p) => {
    if (!p) return;
    if (p.startsWith("**")) out.push({ text: p.slice(2, -2), bold: true });
    else out.push({ text: p, bold: false });
  });
  return out;
}
const plain = (s) => runs(s).map((r) => r.text).join("");
// rough Calibri fit check (avg glyph ≈ 0.48 em)
const WARN = [];
function estLines(text, wIn, pt, em = 0.5) {
  // greedy word wrap with an average Calibri glyph of `em` (long tokens break)
  const cap = (wIn * 72) / (pt * em);
  return String(text).split("\n").reduce((n, l) => {
    let lines = 1, cur = 0;
    for (const word of l.split(/(?<=[\s\-\/–—])/)) {
      let len = word.length;
      if (cur + len > cap && cur > 0) { lines += 1; cur = 0; }
      while (len > cap) { lines += 1; len -= cap; }
      cur += len;
    }
    return n + lines;
  }, 0);
}
function checkFit(label, texts, wIn, hIn, pt, spacing = 1.2, extraPerItem = 0) {
  const lines = texts.reduce((n, t) => n + estLines(t, wIn, pt), 0);
  const need = (lines * pt * spacing) / 72 + texts.length * extraPerItem;
  if (need > hIn) WARN.push(`${label}: needs ~${need.toFixed(2)}in, has ${hIn.toFixed(2)}in`);
}

// ---------- content ----------
const S = []; // slides
let step = 0;
const section = (title, sub) => S.push({ kind: "section", title, sub });
const slide = (o) => { step += 1; S.push({ kind: "content", step, ...o }); };

S.push({
  kind: "title",
  title: "Felodipine ER: everything we did, step by step",
  sub: "MBBE drug 2 · PK-Sim PBPK build, release-rate safe space, fed and dissolution gates, PBPK vs PopPK",
  meta: "PK-Sim, Julia BE simulation and Pumas PopPK",
});

slide({
  title: "The whole project on one page",
  regions: [{ blocks: [
    { type: "flow", steps: [
      { h: "Set-up", t: "CVs, shipped model checks, particle size dropped" },
      { h: "Refit", t: "Weibull refit, disposition check, colonic CYP3A4 = 50%" },
      { h: "Rebuild", t: "Rebuild after lost work, population, coarse 13-run cross" },
      { h: "Safe space", t: "27-run safe space, fed + dissolution gates, PopPK comparison" },
      { h: "Write-up", t: "Summary doc, Word export and this deck" },
    ], h: 1.55 },
    { type: "stats", items: [
      { v: "27", l: "PK-Sim population runs (100 subjects each)" },
      { v: "4.37–5.23 h", l: "ABE safe T50 at reference shape (×0.91–1.09)" },
      { v: "0.79", l: "Model Cmax vs observed (documented bias)" },
      { v: "−7% / −15%", l: "PopPK vs PBPK GMR boundaries" },
    ], h: 1.35 },
    { type: "callout", text: "Bottom line: the fasted T50 × shape safe space is narrow and a diagonal band. Cmax decides BE. PopPK agrees with PBPK near the reference but diverges off-shape.", h: 1.0 },
  ] }],
});

section("1 · Set-up and the model", "Why felodipine, what the OSP model is, and the gates we set before any sweep");

slide({
  title: "Why felodipine was chosen as MBBE drug 2",
  regions: [
    { w: 0.45, blocks: [
      { type: "bullets", items: [
        "**Only OSP library compound that meets all three criteria:** free validated PBPK model, BCS II, and a real marketed ER product (Plendil ER)",
        "**Plendil ER has an FDA product-specific BE guidance**, so the release-rate CQA carries over from metformin",
        "**Goal:** repeat the 2-D T50 × shape safe space on a second, mechanistically opposite drug",
        "**Original extra goal:** particle size and solubility CQAs (dropped, step 6)",
      ] },
    ] },
    { w: 0.55, blocks: [
      { type: "table", fontSize: 13, colW: [2.2, 2.2, 2.3], head: ["", "Metformin XR", "Felodipine ER"], rows: [
        ["BCS class", "III (high sol, low perm)", "II (low sol, high perm)"],
        ["Rate-limiting step", "Permeability, upper-GI window", "In-vivo dissolution"],
        ["Clearance", "Renal, transporters", "Hepatic CYP3A4, high first pass"],
        ["Oral bioavailability", "~50–60%", "~15%"],
        ["Within-subject Cmax CV", "~12%", "~36%"],
        ["Food effect", "Modest", "Large and erratic"],
      ] },
    ] },
  ],
});

slide({
  title: "The OSP felodipine whole-body PBPK model",
  regions: [
    { w: 0.42, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Source:** github.com/Open-Systems-Pharmacology/Felodipine-Model (Bergmann et al., Pharmaceutics 2022;14:1474)",
        "**Built as a CYP3A4 victim drug for DDI work**, not for formulation work; validated at OSP's 2-fold level",
        "**Fitted to 13 studies** (IV, solution, IR, ER) plus an itraconazole DDI",
        "**Structure:** PK-Sim whole body; GI tract as a chain of segments; drug dissolves, crosses the gut wall, meets gut-wall CYP3A4, then the liver",
        "**Our risk:** a DDI model can get AUC right but the absorption shape only roughly",
      ] },
    ] },
    { w: 0.58, blocks: [
      { type: "table", fontSize: 12, colW: [2.15, 2.75, 2.15], head: ["Setting", "Value", "Origin"], rows: [
        ["logP", "4.51", "Fitted (OSP)"],
        ["Fraction unbound", "0.0036", "Ushimura 2010"],
        ["Solubility FaSSIF pH 6.5", "12 µg/mL", "Söderlind 2010"],
        ["CYP3A4 Michaelis–Menten", "Km 2.81 µM; kcat 205 /min", "Walsky 2004; kcat fitted"],
        ["Unspecific hepatic CL", "12.8 /min", "Fitted"],
        ["Gut-wall P (\"P mucosa\")", "0.0892 cm/min", "Fitted, sim-level override"],
        ["Partitioning / permeability", "Rodgers & Rowland / PK-Sim std", "PK-Sim methods"],
        ["ER tablet Weibull", "T50 286.95 min, lag 18.68, shape 0.764", "Fitted"],
        ["IR tablet Weibull", "T50 46.5 min, lag 0, shape 0.888", "Fitted"],
        ["Individual", "European male 30 y, CYP3A4 t½ 36 h", "OSP"],
      ] },
    ] },
  ],
});

slide({
  title: "Four acceptance gates were set before any CQA sweep",
  regions: [{ blocks: [
    { type: "table", fontSize: 14, colW: [0.6, 4.6, 3.0, 3.93], head: ["#", "Gate (scope §3.9)", "Final status", "How it was closed"], rows: [
      ["1", "Oral IR / solution profile reproduced", "Partly", "Solution AUC 1.03; IR AUC overpredicted 1.58× (one dataset)"],
      ["2", "ER fasted Cmax and AUC within 1.25-fold", "Accepted with bias", "AUC 1.04; Cmax 0.79, disposition-driven, documented"],
      ["3", "Fed arm shows direction and size of food effect", "Failed → sensitivity", "Model fed/fasted Cmax 0.96 vs 1.6–2.8 observed"],
      ["4", "Reference Weibull traceable to real dissolution", "Closed, approximate", "Abrahamsson 1994 tablets A/B/C within 1.25-fold"],
    ] },
    { type: "note", text: "Results from the sweep are therefore reported as a fasted claim with a documented Cmax bias and fed as a within-subject-variability sensitivity." },
  ] }],
});

section("2 · First model checks", "Variability inputs, the shipped model against data, and dropping the particle-size CQA");

slide({
  title: "Within-subject CVs sourced from the MHRA public assessment report",
  regions: [
    { w: 0.6, blocks: [
      { type: "table", fontSize: 14, colW: [3.2, 1.6, 2.35], head: ["Study (Delofine XL vs Plendil)", "Cmax CV", "AUC CV"], rows: [
        ["10 mg fasting, 2×2", "35.7%", "24.3% (AUC0–∞)"],
        ["10 mg fed, 2×2", "36.5%", "22.2% (AUC0–∞)"],
        ["10 mg steady state, replicate", "21.6%", "14.9% (AUCτ)"],
        ["2.5 mg fasting, 4-period replicate", "30.8% (CVwR ≈ 34%)", "28.9%"],
      ] },
      { type: "note", text: "Source: UK MHRA PAR, Felodipine/Delofine XL (PL 20117/0397-0399); healthy males." },
    ] },
    { w: 0.4, blocks: [
      { type: "stats", cols: 1, items: [
        { v: "n = 54", l: "ABE 2×2 fasted (Cmax CV 35.7%, GMR 0.95, 80% power)" },
        { v: "n = 56", l: "Fed design (Cmax CV 36.5%)" },
        { v: "n = 24", l: "RSABE / ABEL full replicate, CVwR 34%" },
      ] },
    ] },
  ],
  footnote: "Decision D1: ABE 2-way crossover is primary (FDA guidance); RSABE on Cmax is a sensitivity. Cmax is ~3× metformin's 12% CV; n = 48 is not reused.",
});

slide({
  title: "The shipped model underpredicts ER Cmax by about 22%",
  regions: [
    { w: 0.5, blocks: [
      { type: "table", fontSize: 14, colW: [3.0, 1.4, 1.6], head: ["ER 10 mg study (training)", "Cmax p/o", "AUC p/o"], rows: [
        ["Edgar 1987 (n=12, 12 h)", "0.73", "0.87"],
        ["Bailey 1996 (n=12, 24 h)", "0.63", "1.02"],
        ["Gelal 2005 (n=10, serum)", "1.05", "1.21"],
        ["Bailey 2003 (n=8, 8 h)", "0.76", "1.04"],
        ["**Geometric mean**", "**0.78**", "**1.03**"],
      ] },
    ] },
    { w: 0.5, blocks: [
      { type: "bullets", items: [
        "**AUC passes our 1.25-fold criterion; Cmax fails.** OSP's \"validated\" means within 2-fold, which is looser",
        "**IR check (Edgar 1987, n=12):** Cmax 1.32×, AUC0–12 1.58× — overprediction from 2 h on, so likely disposition, not formulation",
        "**Snapshot formulations confirmed:** ER T50 286.95 min / lag 18.68 / shape 0.764; IR 46.5 / 0 / 0.888; solubility 12 µg/mL",
      ] },
    ] },
  ],
});

slide({
  title: "No particle radius fits the IR data better than the shipped Weibull",
  regions: [
    { w: 0.58, blocks: [
      { type: "table", fontSize: 12, colW: [1.6, 1.2, 1.2, 1.3, 1.15, 1.0], head: ["IR formulation", "Cmax p/o", "tmax (obs 1 h)", "AUC0–12 p/o", "Peaks", "log-RMSE"], rows: [
        ["**Weibull (shipped)**", "1.32", "1.2 h", "1.58", "1", "**0.58**"],
        ["0.5 µm", "2.69", "0.65 h", "1.61", "1", "0.96"],
        ["1 µm", "2.29", "0.85 h", "1.70", "1", "0.85"],
        ["2 µm", "1.68", "1.05 h", "2.00", "1", "0.81"],
        ["3 µm", "1.31", "1.2 h", "2.24", "2", "0.89"],
        ["5 µm", "1.23", "5.75 h", "2.50", "2", "1.02"],
        ["10 µm", "1.11", "7.25 h", "2.42", "1, late", "1.12"],
        ["25 µm", "0.74", "9.5 h", "1.57", "1, late", "1.19"],
        ["50 µm", "0.47", "11.5 h", "0.96", "1, late", "1.42"],
      ] },
    ] },
    { w: 0.42, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Why a double peak:** to lower Cmax, dissolution must slow; undissolved drug reaches the colon, dissolves late and escapes gut CYP3A4, so AUC rises",
        "**Structural trade-off:** no radius avoids both \"Cmax too high\" and \"second peak\"",
        "**PK-Sim 12.2 limit:** Parameter Identification cannot fit particle radius; the \"radius\" entries it offers are physiology — never fit those",
      ] },
    ] },
  ],
});

slide({
  title: "Particle size and solubility were dropped as felodipine CQAs",
  regions: [
    { w: 0.42, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**No public Plendil API particle size** — the patent (US 4,803,081) says only \"micronized\"",
        "**Plendil ER is surfactant-solubilised:** polyoxyl 40 hydrogenated castor oil in an HPMC gel matrix, so drug leaves in solution",
        "**Rat data run opposite to the model:** smaller particles gave higher bioavailability (25 µm 173% vs 200 µm 100%)",
      ] },
    ] },
    { w: 0.58, blocks: [
      { type: "table", fontSize: 12, colW: [0.55, 3.3, 3.2], head: ["#", "Decision", "Why"], rows: [
        ["D1", "MHRA CVs; ABE primary, RSABE sensitivity", "Only public originator ISCV source"],
        ["D2", "No individual-subject exports yet", "Wait for locked reference and population"],
        ["D3", "Drop particle size and solubility", "Structural trade-off; same lever (S/r); not a real-product CQA"],
        ["D4", "Keep T50 × shape, high variability, food", "What felodipine does best, with an FDA guidance"],
        ["D5", "Particle-radius CQA needs another BCS II drug", "Not yet screened"],
        ["D6", "Next: refit the ER Weibull", "ER Cmax 0.78; Cmax decides BE"],
      ] },
    ] },
  ],
});

slide({
  title: "Early export problems led to a checklist used for every export",
  regions: [
    { w: 0.5, blocks: [
      { type: "cards", title: "What went wrong", tone: "accent4", items: [
        "The rad10 export was really the reference sim",
        "rad5 duplicated rad25",
        "A re-copied download \"(1)\" was taken for new results",
        "Two batches ran to 24 h instead of 96 h",
      ] },
    ] },
    { w: 0.5, blocks: [
      { type: "cards", title: "Checklist from then on", tone: "accent5", items: [
        "Configure → Administration: formulation matches the sim name",
        "Output intervals end at 96 h and the chart reaches 96",
        "Unique filename with radius/T50 and \"96h\"",
        "All observed datasets on the chart; I check checksums, end time and parameter-vs-shape",
      ] },
    ] },
  ],
});

slide({
  title: "Data audit: 34 datasets, no grapefruit arms, a train/test split",
  regions: [
    { w: 0.55, blocks: [
      { type: "table", fontSize: 13, colW: [1.6, 5.5], head: ["Group", "Datasets in the OSP snapshot"], rows: [
        ["ER 10 mg train (4)", "Edgar 1987, Bailey 1996, Bailey 2003, Gelal 2005 (serum)"],
        ["ER 10 mg test (6)", "Dresser 2002/2017, Bailey 1998/2000, Lundahl 1997, Madsen 1996"],
        ["ER, other", "Multiple dose: Aberg 1997, Bailey 1995, Lundahl 1998, Blychert 1990; 5/2.5 mg: Goosen 2004, Dresser 2000, Jalava 1997; Guo 2007 (unattached)"],
        ["IR", "Edgar 1987, Smith 1987, Edgar 1992, Bailey 1993, Bratel 1989 (fed), Blychert 1990"],
        ["Solution / IV", "Edgar 1987, Blychert 1990 / Lundahl 1997"],
      ] },
    ] },
    { w: 0.45, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Bailey 1996/2003 and Gelal 2005 are control arms** (observed Cmax 6.1–7.3 nM; an inhibitor arm would be 2–3× higher)",
        "**Correction:** the 0.78 Cmax was on the 4 training studies; the 6 test studies sit in a separate sim; both shipped sims end at 24 h",
        "**Blychert 1990** has ER \"market vs trial\" tablets — a real BE-type comparison for later",
        "**No fed ER data** in the snapshot; Bratel 1989 (IR 5 mg) is the only fed set",
      ] },
    ] },
  ],
});

section("3 · Refit and reference lock", "Weibull refit, disposition diagnosis, colonic CYP3A4, and the locked reference");

slide({
  title: "Against all 10 fasted ER studies the profile is too flat",
  regions: [
    { w: 0.52, blocks: [
      { type: "table", fontSize: 12, colW: [0.8, 2.4, 1.6, 1.8], head: ["Set", "Study", "Cmax p/o", "AUC0–tlast p/o"], rows: [
        ["train", "Edgar 1987", "0.73", "0.87"], ["train", "Bailey 1996", "0.63", "1.07"],
        ["train", "Bailey 2003", "0.76", "1.04"], ["train", "Gelal 2005", "1.05", "1.23"],
        ["test", "Lundahl 1997", "0.83", "1.16"], ["test", "Dresser 2017", "0.76", "0.90"],
        ["test", "Dresser 2002", "1.00", "1.32"], ["test", "Bailey 2000", "0.89", "1.23"],
        ["test", "Bailey 1998", "0.49", "0.81"], ["test", "Madsen 1996", "0.94", "1.43"],
        ["**GM**", "**train / test / all**", "**0.78 / 0.80 / 0.79**", "**1.04 / 1.12 / 1.09**"],
      ] },
    ] },
    { w: 0.48, blocks: [
      { type: "chart", kind: "bar", h: 2.75, title: "Pooled pred/obs by time window", labels: ["0–2 h", "2–6 h", "6–12 h", "12–24 h"],
        series: [{ name: "pred/obs", values: [1.23, 0.90, 1.22, 3.18] }], yMin: 0, yMax: 3.5, valueFmt: "0.00" },
      { type: "bullets", fontSize: 13, items: [
        "**Too much drug early, too little at the peak, too much in the tail** — the signature of Weibull shape < 1",
        "Predicted Cmax 4.60 nM at 3.0 h vs observed 4.4–9.4 nM at 2.5–5 h",
      ] },
    ] },
  ],
});

slide({
  title: "The Weibull refit worked as a fit but did not fix Cmax",
  regions: [
    { w: 0.4, blocks: [
      { type: "stats", cols: 1, items: [
        { v: "T50 296 ± 52 min", l: "95% CI ≈ 244–348; shipped 287 inside" },
        { v: "Shape 0.89 ± 0.08", l: "≈ 0.81–0.97; shipped 0.764 outside" },
        { v: "r = −0.53", l: "T50–shape correlation: not confounded" },
      ] },
      { type: "note", text: "PI on the 4 training studies, log scaling, Levenberg–Marquardt, 10 random starts, lag fixed; all starts reached total error 0.95." },
    ] },
    { w: 0.6, blocks: [
      { type: "table", fontSize: 13, colW: [1.7, 1.05, 1.05, 0.85, 0.85, 0.85, 0.9], head: ["", "Cmax GM", "AUC GM", "0–2 h", "2–6 h", "6–12 h", "12–24 h"], rows: [
        ["train, shipped", "0.78", "1.04", "1.31", "0.86", "1.19", "1.83"],
        ["train, refit", "0.77", "1.07", "1.02", "0.87", "1.31", "1.82"],
        ["test, shipped", "0.80", "1.12", "1.18", "0.92", "1.24", "4.60"],
        ["test, refit", "0.79", "1.14", "0.91", "0.94", "1.36", "4.74"],
      ] },
      { type: "callout", text: "The refit fixes the early phase, slightly worsens the tail and leaves Cmax at 0.77–0.79. Observed profiles peak higher and fall faster at the same AUC — a release change cannot do that, so the cause is disposition." },
    ] },
  ],
});

slide({
  title: "IV and solution: disposition is off, oral AUC right by cancellation",
  regions: [
    { w: 0.5, blocks: [
      { type: "table", fontSize: 13, colW: [2.2, 1.3, 1.6, 1.0], head: ["IV 1.5 mg, Lundahl 1997", "Model", "Observed", "Ratio"], rows: [
        ["AUC∞ (µmol·min/L)", "3.63", "4.59", "0.79"],
        ["Cmax end of infusion (nM)", "22.3", "30.0", "0.74"],
        ["CL (mL/min/kg)", "14.7", "11.6", "**1.27**"],
        ["Vss (mL/kg)", "4985", "3090", "**1.61**"],
        ["t½ (h)", "8.1", "5.7", "**1.42**"],
      ] },
      { type: "note", text: "Oral solution 10 mg (Edgar 1987): AUC∞ 1.03, Cmax 1.24. Implied F: model 0.20 vs observed 0.15 (literature ~15%)." },
    ] },
    { w: 0.5, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**CL too high (×1.27) and F too high (×1.3) cancel** in oral AUC",
        "**Excess Vss and t½ flatten slow-input profiles:** low peak, high tail — most of the ER Cmax 0.79",
        "**ER tail excess (1.8–4.7) is bigger than IV or solution tails**, so a second, ER-specific cause is likely: colonic absorption without gut first pass",
      ] },
      { type: "cards", title: "Options considered", tone: "accent1", items: [
        "A — accept the documented Cmax bias, use only T/R ratios",
        "B — recalibrate disposition (risks OSP's DDI validation)",
        "**C — colon-specific lever, then T50 refit (chosen)**",
      ] },
    ] },
  ],
});

slide({
  title: "Colonic CYP3A4 at 50% of the small intestine fixes the tail, not Cmax",
  regions: [
    { w: 0.56, blocks: [
      { type: "table", fontSize: 13, colW: [1.9, 1.0, 1.15, 1.15, 1.9], head: ["Colonic CYP3A4", "F", "Cmax GM", "AUC GM", "12–24 h train / test"], rows: [
        ["0 (shipped)", "~0.20", "0.789", "1.089", "1.83 / 4.60"],
        ["10% of SI (0.0073)", "0.200", "0.789", "1.068", "1.49 / 3.93"],
        ["25% of SI (0.018)", "0.186", "0.788", "1.054", "1.27 / 3.48"],
        ["**50% of SI (0.036)**", "**0.176**", "**0.787**", "**1.042**", "**1.12 / 3.14**"],
      ] },
      { type: "note", text: "Shipped model: CYP3A4 in liver and small-intestine segments (0.073) but none in the large intestine, while AADAC and P-gp do have colon entries." },
    ] },
    { w: 0.44, blocks: [
      { type: "cards", title: "Reference locked (decision)", tone: "accent1", items: [
        "Shipped ER Weibull: T50 286.95 min, lag 18.68 min, shape 0.764",
        "Colonic CYP3A4 = 0.036 in all large-intestine segments",
        "Second refit (step 5c) skipped",
        "Cmax bias ≈ 0.79 documented as disposition-driven",
      ] },
      { type: "note", text: "Logged as a deviation from scope §3.9: gates #2–#4 not yet met, so sweep results were provisional until they closed." },
    ] },
  ],
});

slide({
  title: "Rebuild did not match, and the work was not saved",
  regions: [
    { w: 0.5, blocks: [
      { type: "table", fontSize: 14, colW: [2.0, 2.0, 2.0], head: ["", "Colon-50 clone", "New individual"], rows: [
        ["Cmax", "4.59 nM @ 3.0 h", "4.72 nM @ 2.75 h"],
        ["F", "0.176", "0.159"],
        ["Fraction absorbed", "0.981", "0.959"],
        ["CL (mL/min/kg)", "14.27", "14.35"],
        ["0–2 h p/o", "1.31", "1.47"],
      ] },
    ] },
    { w: 0.5, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Stage 6a:** new profile `CYP3A4 36 h colon50` and individual `European_colonCYP3A4_50`, sim switched via Configure — check failed",
        "**Suspected cause:** OSP stores 27 simulation-level overrides that Configure may drop",
        "**Later finding:** none of this PK-Sim work was saved (clones, profile, individual); only the xlsx exports survived",
        "**Rule since:** always end a PK-Sim stage with File → Save",
      ] },
    ] },
  ],
});

section("4 · Rebuild and population", "The override hunt, the virtual population, BE code, and the coarse cross");

slide({
  title: "Configure drops OSP's simulation overrides; Clone keeps them",
  regions: [
    { w: 0.55, blocks: [
      { type: "table", fontSize: 12, colW: [1.6, 4.0, 1.5], head: ["Step", "What we found", "Result"], rows: [
        ["9.8", "Original individual checked: LI CYP3A4 = 0, SI 0.07", "Clean"],
        ["9.9", "First snapshot sent was the older particle-size project", "Wrong file"],
        ["9.10", "22 P(int↔cell) = OSP fitted \"P mucosa\" 0.0892 cm/min", "Not cosmetic"],
        ["9.11", "Rebuild without re-entry = earlier failure to the digit", "Overrides absent"],
        ["9.12", "Screenshots: P mucosa survived Configure (User Defined)", "Refuted for 22"],
        ["9.13", "Snapshot diff: only Tablet time delay factor missing", "Cause found"],
        ["9.14", "v2 with delay factor 0.0695415108 re-entered", "Match to 3.3e-5"],
      ] },
    ] },
    { w: 0.45, blocks: [
      { type: "cards", title: "The 27 OSP overrides (reference/ER_tablet_10mg_sim_overrides.csv)", tone: "accent1", items: [
        "22 × gut-wall P (interstitial↔intracellular) = 0.0891674187 cm/min",
        "Tablet time delay factor = 0.0695415108",
        "Water volume 3.5 mL/kg; Undefined Liver × 3",
      ] },
      { type: "callout", text: "Configure → individual drops the delay factor. Configure → population also drops the 22 P mucosa values. Clone keeps everything." },
    ] },
  ],
  footnote: "Reference restored in Felodipine_ER_colon50.pksim5: Cmax 4.588 nM @ 3.0 h, AUC0–96 4.383 µM·min, F 0.1765, all-10 Cmax 0.787 / AUC 1.042.",
});

slide({
  title: "One 100-subject population, reused for every formulation",
  regions: [
    { w: 0.42, blocks: [
      { type: "table", fontSize: 13, colW: [1.7, 3.6], head: ["Pop_felo_100", "Value"], rows: [
        ["Base individual", "European_colonCYP3A4_50"],
        ["Size / sex", "100 subjects, 50 F / 50 M"],
        ["Age", "18.2–54.7 y (requested 18–55)"],
        ["Weight / BMI", "41–98 kg / 19–30 kg/m²"],
        ["Dose", "10 mg ER, single dose, fasted"],
        ["Run length", "0–96 h"],
        ["Overrides", "23 re-entered by hand in Pop_ref"],
      ] },
      { type: "note", text: "100, not 24: trials need ~54 subjects and resampling from 24 physiologies duplicates heavily. Every sweep run is a clone of Pop_ref, so test and reference are paired." },
    ] },
    { w: 0.58, blocks: [
      { type: "table", fontSize: 13, colW: [1.7, 1.15, 1.05, 1.0, 1.35, 0.85], head: ["Pop_ref", "Individual", "Median", "GM", "Range", "BSV CV"], rows: [
        ["Cmax (nM)", "4.59", "5.26", "5.16", "1.65–10.98", "35%"],
        ["tmax (h)", "3.0", "2.75", "2.71", "1.9–3.75", "15%"],
        ["AUC0–96 (µM·min)", "4.38", "5.18", "5.16", "1.39–12.98", "39%"],
      ] },
      { type: "bullets", fontSize: 14, items: [
        "**~1.15× above the individual, uniformly over time** — fits the demographics (half lighter women, fixed dose); not the missing-delay signature (that is early up, late down)",
        "**Accepted as a correct build**",
      ] },
    ] },
  ],
});

slide({
  title: "be_utils.jl was refactored so one BE engine serves both drugs",
  regions: [
    { w: 0.55, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Metformin API unchanged:** all 5 metformin output CSVs byte-identical before and after",
        "**Added:** DRUG_CONFIG, wsv_cv(drug, arm, endpoint), parse_formulation (T50, shape), drug-aware loader (venous blood rows only, errors on duplicates)",
        "**Added:** power_abe_2x2 / sample_size_abe_2x2, full-replicate FDA RSABE (sWR ≥ 0.294) and EMA ABEL (cap at CVwR 50%)",
        "**21 tests pass** (programs/julia/test_be_utils.jl)",
      ] },
    ] },
    { w: 0.45, blocks: [
      { type: "table", fontSize: 13, colW: [3.6, 1.8], head: ["Test", "Result"], rows: [
        ["n at CV 0.20 / 0.25 / 0.30", "20 / 28 / 40 (PowerTOST)"],
        ["Felodipine fasted Cmax CV 0.357", "n = 54"],
        ["Felodipine fed Cmax CV 0.365", "n = 56"],
        ["Loader reproduces T50 2.39 h GMR", "1.488"],
        ["RSABE self-comparison, n = 24", "97% pass"],
        ["True GMR 1.25: ABE / RSABE pass", "4.6% / 20.8%"],
      ] },
      { type: "note", text: "The 20.8% at GMR 1.25 is the known RSABE type I inflation — worth reporting." },
    ] },
  ],
});

slide({
  title: "Coarse cross (13 runs): Cmax moves with T50 and is U-shaped in shape",
  regions: [
    { w: 0.5, blocks: [
      { type: "chart", kind: "line", h: 3.6, title: "Cmax GMR vs T50 (shape 0.764)", labels: ["2.39", "3.35", "4.07", "4.78", "5.74", "7.17", "9.56"], xTitle: "T50 (h)",
        series: [{ name: "Cmax GMR", values: [1.488, 1.229, 1.099, 1.0, 0.898, 0.786, 0.659] },
                 { name: "1.25", values: [1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25], dash: true },
                 { name: "0.80", values: [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8], dash: true }], yMin: 0.6, yMax: 1.6, valueFmt: "0.00" },
    ] },
    { w: 0.5, blocks: [
      { type: "chart", kind: "line", h: 3.6, title: "Cmax GMR vs shape (T50 4.78 h)", labels: ["0.5", "0.6", "0.764", "0.9", "1.0", "1.2", "1.5"], xTitle: "Weibull shape",
        series: [{ name: "Cmax GMR", values: [1.229, 1.094, 1.0, 1.010, 1.037, 1.119, 1.269] },
                 { name: "1.25", values: [1.25, 1.25, 1.25, 1.25, 1.25, 1.25, 1.25], dash: true },
                 { name: "0.80", values: [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8], dash: true }], yMin: 0.6, yMax: 1.6, valueFmt: "0.00" },
    ] },
  ],
  footnote: "Joint ABE power ≥ 0.80 (n = 54) only inside T50 ×0.85–1.2; shape tolerates 0.764→1.0 but not 0.6 or 1.2. AUC passes almost everywhere. Every T50/shape was verified from the lumen fraction dissolved.",
});

section("5 · The safe space", "2-D fill, edge runs, and the 27-formulation result");

slide({
  title: "27 formulations were simulated in three rounds",
  regions: [
    { w: 0.62, blocks: [{ type: "image", src: "slides/full-summary/fig_design_nodates.png" }] },
    { w: 0.38, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Coarse cross (12 + ref):** T50 axis and shape axis",
        "**2-D fill (10):** T50 4.07 and 5.74 h × shapes 0.6–1.2, plus 4.42 and 5.26 h at 0.764",
        "**Band extension (4):** 6.31 h × 0.6/1.0/1.2 and 6.88 h × 1.2",
        "**Upload checks caught:** a missing PK-Analyses file (re-exported) and \"5.7 h\" files that were really 344 min (renamed 5.74 h)",
      ] },
    ] },
  ],
});

const RUNS = [
  ["2.39 (143)", "0.764", "1.488", "1.040", "0.00", "0.99", "0.00", "0.00", "0.00"],
  ["3.35 (201)", "0.764", "1.229", "1.026", "0.10", "1.00", "0.08", "0.29", "0.20"],
  ["4.07 (244)", "0.6", "1.189", "0.970", "0.19", "0.99", "0.18", "0.45", "0.33"],
  ["4.07 (244)", "0.764", "1.099", "1.013", "0.62", "1.00", "0.59", "0.81", "0.73"],
  ["4.07 (244)", "0.9", "1.116", "1.032", "0.52", "0.99", "0.49", "0.75", "0.66"],
  ["4.07 (244)", "1.0", "1.149", "1.040", "0.35", "0.99", "0.32", "0.62", "0.52"],
  ["4.07 (244)", "1.2", "1.243", "1.050", "0.07", "0.98", "0.06", "0.25", "0.14"],
  ["4.42 (265)", "0.764", "1.047", "1.007", "0.84", "1.00", "0.83", "0.94", "0.91"],
  ["4.78 (287)", "0.5", "1.229", "0.904", "0.10", "0.80", "0.06", "0.29", "0.19"],
  ["4.78 (287)", "0.6", "1.094", "0.951", "0.64", "0.97", "0.61", "0.82", "0.74"],
  ["**4.78 (287) ref**", "**0.764**", "**1.000**", "**1.000**", "0.90", "1.00", "**0.90**", "0.97", "0.95"],
  ["4.78 (287)", "0.9", "1.010", "1.022", "0.89", "1.00", "0.89", "0.97", "0.94"],
  ["4.78 (287)", "1.0", "1.037", "1.032", "0.86", "0.99", "0.85", "0.95", "0.92"],
  ["4.78 (287)", "1.2", "1.119", "1.045", "0.50", "0.98", "0.46", "0.74", "0.65"],
  ["4.78 (287)", "1.5", "1.269", "1.055", "0.03", "0.97", "0.03", "0.17", "0.10"],
  ["5.26 (316)", "0.764", "0.946", "0.991", "0.79", "1.00", "0.79", "0.91", "0.87"],
  ["5.74 (344)", "0.6", "0.996", "0.927", "0.90", "0.91", "0.82", "0.97", "0.94"],
  ["5.74 (344)", "0.764", "0.898", "0.981", "0.52", "1.00", "0.54", "0.76", "0.67"],
  ["5.74 (344)", "0.9", "0.901", "1.008", "0.53", "1.00", "0.55", "0.76", "0.68"],
  ["5.74 (344)", "1.0", "0.922", "1.020", "0.67", "1.00", "0.67", "0.85", "0.79"],
  ["5.74 (344)", "1.2", "0.990", "1.036", "0.90", "0.99", "0.88", "0.97", "0.94"],
  ["6.31 (379)", "0.6", "0.947", "0.913", "0.79", "0.85", "0.69", "0.92", "0.88"],
  ["6.31 (379)", "1.0", "0.865", "1.012", "0.30", "1.00", "0.31", "0.58", "0.44"],
  ["6.31 (379)", "1.2", "0.927", "1.030", "0.70", "0.99", "0.70", "0.88", "0.81"],
  ["6.88 (413)", "1.2", "0.873", "1.023", "0.35", "1.00", "0.36", "0.62", "0.49"],
  ["7.17 (430)", "0.764", "0.786", "0.952", "0.03", "0.97", "0.03", "0.14", "0.09"],
  ["9.56 (574)", "0.764", "0.659", "0.904", "0.00", "0.81", "0.00", "0.00", "0.00"],
];
const RUN_HEAD = ["T50 h (min)", "Shape", "BA Cmax GMR", "BA AUC0–t GMR", "ABE Cmax", "ABE AUC", "ABE joint", "RSABE", "ABEL"];
const RUN_COLW = [1.85, 0.95, 1.45, 1.6, 1.25, 1.2, 1.3, 1.15, 1.1];
const passRow = (r) => parseFloat(plain(r[6])) >= 0.8;
slide({
  title: "Matrix sweep results, runs 1–14: relative BA and power",
  regions: [{ blocks: [
    { type: "table", fontSize: 12, colW: RUN_COLW, head: RUN_HEAD, rows: RUNS.slice(0, 14), highlight: passRow },
  ] }],
  footnote: "BA = paired test/reference geometric-mean ratio in 100 subjects. Power = share of 1000 simulated trials passing: ABE 2×2 n = 54 (CV Cmax 35.7%, AUC0–t 25.3%); RSABE / ABEL on Cmax, n = 24, CVwR 34%. Shaded: joint power ≥ 0.80.",
});
slide({
  title: "Matrix sweep results, runs 15–27: relative BA and power",
  regions: [{ blocks: [
    { type: "table", fontSize: 12, colW: RUN_COLW, head: RUN_HEAD, rows: RUNS.slice(14), highlight: passRow },
  ] }],
  footnote: "AUC0–t power is ≥ 0.85 everywhere except 4.78 h / shape 0.5 and 9.56 h (0.80–0.81): Cmax is the binding metric. Scripts felo_02_be_sweep.jl → output_sweep_2d/.",
});

slide({
  title: "The safe space is a diagonal band, not a box",
  regions: [
    { w: 0.6, blocks: [{ type: "image", src: "output_sweep_2d/fig_safespace_2d_joint.png" }] },
    { w: 0.4, blocks: [
      { type: "table", fontSize: 12, colW: [0.85, 1.35, 1.35, 1.25], head: ["Shape", "ABE n=54", "RSABE n=24", "ABEL n=24"], rows: [
        ["0.6", "5.25–5.84", "4.73–≥6.31", "4.90–≥6.31"],
        ["**0.764**", "**4.37–5.23**", "4.07–5.63", "4.19–5.46"],
        ["0.9", "4.47–5.28", "4.16–5.66", "4.29–5.49"],
        ["1.0", "4.67–5.49", "4.36–5.87", "4.49–5.70"],
        ["1.2", "5.24–6.10", "4.90–6.50", "5.04–6.33"],
      ] },
      { type: "bullets", fontSize: 13, items: [
        "T50 range (h) with joint power ≥ 0.80",
        "**At fixed T50, Cmax is lowest near shape 0.8**, so a slower T50 is offset by moving shape either way",
        "**Method:** log GMR interpolated on a Delaunay mesh, isotonic calibration to power (within 0.036 of all 27 runs)",
      ] },
    ] },
  ],
});

section("6 · Fed and dissolution gates", "Gate #3 food effect and gate #4 in-vitro link");

slide({
  title: "Fed literature sets the target: Cmax up 1.6–2.8×, AUC unchanged",
  regions: [
    { w: 0.55, blocks: [
      { type: "table", fontSize: 13, colW: [2.6, 2.0, 2.5], head: ["Source", "Fed / fasted Cmax", "Note"], rows: [
        ["Sandoz CA monograph (Plendil arm)", "×2.77", "AUC∞ ×1.09; tmax 4.5 h both"],
        ["US Plendil label", "≈ +60%", "AUC unchanged"],
        ["Weitschies 2005 (n=6)", "×1.63", "Bursts from proximal stomach"],
        ["Apotex CA monograph", "×3.54", "Different subjects; not used"],
      ] },
      { type: "note", text: "Profiles: Kim 2014 PLoS One figs (digitisable, < 8 h); Abrahamsson 1993 Pharm Res (paywalled). Registered fed BE trials (NCT02311530, NCT01630902) posted no results. Full write-up: reference/fed-arm-literature.md." },
    ] },
    { w: 0.45, blocks: [
      { type: "stats", cols: 1, items: [
        { v: "1.6–2.8", l: "Target Cmax fed/fasted" },
        { v: "1.0–1.1", l: "Target AUC fed/fasted" },
        { v: "0.81 / 0.89", l: "Fasted model vs Sandoz Plendil: Cmax / AUC∞" },
      ] },
    ] },
  ],
});

slide({
  title: "The model cannot produce the food effect, so fed became a sensitivity",
  regions: [
    { w: 0.48, blocks: [
      { type: "stats", cols: 2, items: [
        { v: "0.96", l: "Cmax fed/fasted (target 1.6–2.8)" },
        { v: "1.01", l: "AUC∞ fed/fasted (target ≈ 1.0–1.1)" },
        { v: "2.75 → 5.0 h", l: "Median tmax (observed unchanged at 4.5 h)" },
        { v: "Unchanged", l: "Safe space with fed CV 36.5%, n = 56" },
      ] },
      { type: "note", text: "Pop_fed_ref = clone of Pop_ref + high-fat breakfast event at t = 0; Configure kept the P mucosa values this time. The meal template used was not recorded." },
    ] },
    { w: 0.52, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**Why it fails:** the Weibull release is a fixed function of time, so the meal only slows gastric emptying — Cmax falls slightly",
        "**What really happens:** the tablet sits in the proximal stomach, released drug empties in bursts and more escapes gut-wall extraction (Weitschies 2005, Bergstrand 2009)",
        "**PK-Sim's meal event has neither mechanism**",
        "**Decision (option A):** claim the fasted safe space only; fed = within-subject-variability sensitivity (scope §10)",
      ] },
    ] },
  ],
});

slide({
  title: "Dissolution link: the in-vivo Weibull fails USP Test 1 at 10 h",
  regions: [
    { w: 0.55, blocks: [
      { type: "table", fontSize: 13, colW: [3.3, 1.2, 1.2, 1.4], head: ["USP Felodipine ER Test 1", "2 h", "6 h", "10 h"], rows: [
        ["Limits (App. 2, 50 rpm, pH 6.5 + 1% SLS)", "10–30%", "42–68%", "≥ 75%"],
        ["Reference in-vivo Weibull, read as in vitro", "26.9%", "54.7%", "**69.6% — fails**"],
      ] },
      { type: "bullets", fontSize: 14, items: [
        "**Real Plendil passes USP**, so in vitro and in vivo differ in shape, not just timescale; a 1:1 map is wrong",
        "**Abrahamsson 1994 (Pharm Res):** tablets A/B/C on the same HPMC platform, n = 16, IVIVC r² 0.94–0.97",
        "**Its SLS 50 rpm row = USP Test 1:** reference MDT 9.4 h sits inside the USP-implied 4.5–12.9 h — consistent, but a weak constraint",
        "Vetchy 2007 (no basket) cannot be chained to the correlation",
      ] },
    ] },
    { w: 0.45, blocks: [
      { type: "cards", title: "Proposal that followed", tone: "accent1", items: [
        "Fit Weibulls to the paper's in-vivo dissolution for tablets A, B and C",
        "Run each as a clone of Pop_ref in PK-Sim",
        "Compare Cmax levels and Cmax ratios with the observed mean curves",
        "Use it as an external check of the model's Cmax sensitivity to release rate",
      ] },
    ] },
  ],
});

slide({
  title: "Tablets A/B/C: three real formulations reproduced within 1.25-fold",
  regions: [
    { w: 0.5, blocks: [
      { type: "chart", kind: "bar", h: 3.3, title: "Mean-curve peak (nM)", labels: ["A (102 min, b 1.18)", "B (363 min, b 0.94)", "C (1185 min, b 1.10)"],
        series: [{ name: "Model", values: [11.4, 4.8, 2.1] }, { name: "Observed", values: [9.5, 5.0, 2.3] }], yMin: 0, yMax: 12, valueFmt: "0.0", legend: true },
      { type: "note", text: "Model / observed: A 1.20, B 0.96, C 0.91." },
    ] },
    { w: 0.5, blocks: [
      { type: "table", fontSize: 13, colW: [1.6, 1.3, 1.4, 1.5], head: ["Ratio", "Model", "Observed", "Difference"], rows: [
        ["Cmax B/C", "2.31", "2.17", "+6%"],
        ["Cmax A/B", "2.36", "1.90", "+24%"],
        ["AUC C/B", "0.80", "0.93", "slow tablet under-absorbed"],
      ] },
      { type: "bullets", fontSize: 14, items: [
        "**Slow side right (+6%); fast side overstated ~25%**, so the band's short-T50 edge is conservative",
        "**Gate #4 closed** as \"consistent, approximate, with a documented fast-side conservatism\"",
        "**Caveats:** Weibulls from deconvolution, digitised by eye; observed are young Swedish men, n = 16",
      ] },
    ] },
  ],
});

section("7 · PBPK vs PopPK", "Same pipeline as metformin: fit PopPK to simulated studies and compare BE decisions");

slide({
  title: "PopPK set-up: per-arm fits on rich simulated studies",
  regions: [
    { w: 0.52, blocks: [
      { type: "flow", vertical: true, steps: [
        { h: "01 dataset", t: "PBPK subjects → 18 samples (0.5–96 h), 15% proportional error; 27 arms × 100" },
        { h: "02 fit", t: "Per-arm FOCE, both basin starts, keep higher LL; NCA on 0.05 h IPRED grid" },
        { h: "03 collect", t: "Copy to poppk/output_rich/ in PK-Sim format so felo_02/03 run unchanged" },
        { h: "04–06", t: "Compare + type I (20k trials), diagnostics, design figure" },
      ] },
    ] },
    { w: 0.48, blocks: [
      { type: "bullets", fontSize: 14, items: [
        "**2-cmt + Weibull was unidentifiable** (Q or Vp hit bounds), so 1-cmt + Weibull, lag 0.3113 h, IIV on alpha, CL, Vc",
        "**Two likelihood basins:** release-limited (flip-flop) and disposition-limited (Vc ≈ 4000 L, t½ ≈ 9 h)",
        "**23 arms disposition-limited, 4 release-limited;** σ 0.16–0.24 vs true 0.15",
        "**PopPK release parameters do not map onto PK-Sim T50**, unlike metformin",
        "Run: julia -J pumas.so --project=julia-env",
      ] },
    ] },
  ],
});

slide({
  title: "PopPK agrees near the reference but the band diverges off-shape",
  regions: [
    { w: 0.58, blocks: [{ type: "image", src: "output_sweep_2d/fig_safespace_2d_overlay_PopPK.png" }] },
    { w: 0.42, blocks: [
      { type: "table", fontSize: 12, colW: [0.9, 1.75, 2.15], head: ["Shape", "PopPK (h)", "PBPK (h)"], rows: [
        ["0.6", "3.91–4.79", "5.25–5.84 (no overlap)"],
        ["0.764", "4.32–5.25", "4.37–5.23"],
        ["0.9", "4.67–5.60", "4.47–5.28"],
        ["1.0", "4.77–5.71", "4.67–5.49"],
        ["1.2", "5.10–6.11", "5.24–6.10"],
      ] },
      { type: "bullets", fontSize: 13, items: [
        "**Per subject:** AUC0–t within ±3% in every arm; Cmax 0.87 (shape 0.5) to 1.135; ref 1.089",
        "**Bias does not cancel in T/R:** shape 0.5 Cmax GMR 1.229 PBPK vs 0.984 PopPK",
        "**1-D T50 axis agrees within 0.08 h** (4.38–5.21 vs 4.31–5.21)",
      ] },
    ] },
  ],
});

slide({
  title: "Type I error: both are calibrated, but PopPK moves the GMR boundaries",
  regions: [
    { w: 0.6, blocks: [
      { type: "table", fontSize: 13, colW: [1.6, 1.3, 1.3, 1.25, 1.95], head: ["T50 axis", "PBPK", "PopPK", "Shift", "True GMR there / truth passes"], rows: [
        ["Fast side (GMR 1.25)", "3.25 h", "2.78 h", "−15%", "1.37 / 0.2% (n=54), 0.5% (n=24)"],
        ["Slow side (GMR 0.80)", "6.96 h", "6.45 h", "−7%", "0.84 / 15.8% (n=54), 8.2% (n=24)"],
      ] },
      { type: "bullets", fontSize: 14, items: [
        "**Calibrated at their own boundary:** 4.5–5.1% for Cmax and joint at n = 54; TOST control 4.96%",
        "**Fast side anti-conservative, slow side conservative**",
        "**Shape axis:** PBPK reaches GMR 1.25 at 1.46; PopPK never does",
        "**RSABE / ABEL at the ABE boundary (n=24):** 20–22% / 13–14% for both",
      ] },
    ] },
    { w: 0.4, blocks: [
      { type: "stats", cols: 1, items: [
        { v: "−1.1% / −6.8%", l: "Metformin PopPK vs PBPK boundaries" },
        { v: "−7% / −15%", l: "Felodipine PopPK vs PBPK boundaries" },
      ] },
      { type: "callout", text: "PopPK–PBPK decision agreement needs structural adequacy; for felodipine it fails off-reference even with rich sampling." },
    ] },
  ],
  footnote: "Discussion deck (18 slides): docs/Discussion-deck_PBPK-vs-PopPK_offline.html · scripts felo_04 (compare, type I), felo_05 (diagnostics), felo_06 (design).",
});

section("8 · Wrap-up", "Lessons, where everything lives, and what is still open");

slide({
  title: "PK-Sim gotchas we learned the hard way",
  regions: [
    { w: 0.5, blocks: [
      { type: "cards", title: "Model handling", tone: "accent4", items: [
        "Configure drops OSP sim-level overrides; Clone keeps them — re-enter from the overrides CSV",
        "Overrides show under Parameters → User Defined (+ Find box); neighbourhoods in the detailed tree",
        "Do not \"update from building block\" on colon clones — it resets colonic CYP3A4 to 0",
        "PI cannot fit particle radius; never fit the physiology \"radius\" entries",
      ] },
    ] },
    { w: 0.5, blocks: [
      { type: "cards", title: "Working remotely", tone: "accent1", items: [
        "Always File → Save at the end of each stage (a full day of work was lost once)",
        "Verify every export: checksums, 96 h end, observed data, parameter vs curve",
        "Ask for a screenshot instead of guessing a UI label",
        "Ask for a project snapshot JSON — and check it is the current project",
      ] },
    ] },
  ],
});

slide({
  title: "Where everything is saved",
  regions: [{ blocks: [
    { type: "table", fontSize: 12, colW: [3.4, 8.73], head: ["Item", "Location (Felodipine-BCSII/ unless stated)"], rows: [
      ["PK-Sim project (local installation)", "Felodipine_ER_colon50.pksim5 (baseline Felodipine_original.pksim5 untouched)"],
      ["Project notes, all decisions", "project notes (not included) · scope: planning notes (not included)"],
      ["PK-Sim exports", "pksim-particle-size-pk/ · pksim-releaserate-cqa-er/ (matrix_Sweep_td50_shape/, 16Oct2026sims/, fed-sims/, tablets-1994/)"],
      ["Sweep results", "output_sweep_coarse/ · output_sweep_2d/ · output_sweep_2d_fedcv/ · output_tablets1994/"],
      ["PopPK", "poppk/ (src, scripts, output_rich/, output_sweep_2d_poppk/) · runs in poppk-fits/"],
      ["Analysis scripts", "scripts/felo_02 … felo_06 (run from programs, julia --project=.) · scripts/score.py, xl.py"],
      ["References", "reference/: OSP snapshot JSON, overrides CSV, fed + dissolution literature, digitised/"],
      ["Decks", "Explainer v14: docs/Explainer-deck_PKSim-model-fix_offline.html · Discussion: docs/Discussion-deck_PBPK-vs-PopPK_offline.html"],
      ["Summary doc", "docs/Felodipine-ER-PKSim-Work-Summary.docx · Word copy: Felodipine-ER-PKSim-Work-Summary.docx"],
      ["This deck", "slides/full-summary/Felodipine-ER-Full-Work-Summary.pptx and .html · generator scripts/build_full_deck.js"],
    ] },
  ] }],
});

slide({
  title: "PK-Sim work is finished; these decisions are open",
  regions: [
    { w: 0.5, blocks: [
      { type: "cards", title: "Proposed next steps", tone: "accent1", items: [
        "NCA on the sampled data (model-free comparator)",
        "Joint PopPK model across arms (T50 and shape covariates)",
        "Richer absorption model (transit or two-phase)",
        "Sparse sampling (6 samples)",
        "Felodipine write-up",
      ] },
    ] },
    { w: 0.5, blocks: [
      { type: "cards", title: "Open questions and loose ends", tone: "accent2", items: [
        "Report the PopPK divergence as is, or improve PopPK first?",
        "Literature: colon vs SI CYP3A4 ratio (supports 50%); Blychert 1990 market vs trial",
        "Record the meal template used in Pop_fed_ref",
        "Class B (particle size) CQA needs a different BCS II drug — not screened",
        "Refresh the MBBE archive with the felodipine work when asked",
      ] },
    ] },
  ],
});

// ---------- layout ----------
// Each region (column) stacks its blocks. Text blocks get the largest font bump (0–4 pt)
// whose estimated height still fits; images, charts and stat columns share what is left.
const BASE_FS = { table: 12, bullets: 15, cards: 14 };
const MAX_FS = { table: 16, bullets: 20, cards: 19 };
function fsOf(b) { return b._fs || b.fontSize || BASE_FS[b.type] || 12; }
function estH(b, w, adj) {
  const base = b.fontSize || BASE_FS[b.type];
  const fsz = base ? Math.min(base + adj, Math.max(base, MAX_FS[b.type])) : 0;
  if (b.type === "table") {
    const colW = scaleCols(b.colW, w) || b.head.map(() => w / b.head.length);
    const rowLines = [b.head, ...b.rows].map((r, ri) => Math.max(...r.map((c, i) => estLines(plain(c), colW[i] - 0.14, fsz, ri === 0 || String(c).startsWith("**") ? 0.56 : 0.46))));
    return { fs: fsz, h: 1.04 * rowLines.reduce((a, n) => a + (fsz * 1.2 * n) / 72 + 0.085, 0) };
  }
  if (b.type === "bullets") {
    const lines = b.items.reduce((a, t) => a + estLines(plain(t), w - 0.3, fsz), 0);
    return { fs: fsz, h: (lines * fsz * 1.2) / 72 + b.items.length * 0.1 + 0.05 };
  }
  if (b.type === "cards") {
    const lines = b.items.reduce((a, t) => a + estLines(plain(t), w - 0.75, fsz), 0);
    return { fs: fsz, h: 0.8 + (lines * fsz * 1.2) / 72 + b.items.length * 0.1 + 0.1 };
  }
  if (b.type === "note") return { h: Math.max(0.3, estLines(b.text, w, 11) * 0.19 + 0.05) };
  if (b.type === "callout") return { h: estLines(plain(b.text), w - 0.4, 15) * 0.26 + 0.3 };
  if (b.type === "stats" && (b.cols || b.items.length) >= b.items.length) return { h: b.h || 1.3 };
  if (b.h) return { h: b.h };
  return null; // flexible: image, chart, stacked stats, vertical flow
}
const MIN_FLEX = { image: 2.6, chart: 2.6, stats: 2.4, flow: 3.0 };
function layoutSlide(s) {
  if (s._placed) return s._placed;
  const top = CONTENT.y, avail = CONTENT.h - (s.footnote ? 0.6 : 0);
  const regs = s.regions || [];
  const totalGap = GAP * (regs.length - 1);
  let x = CONTENT.x;
  const placed = [];
  regs.forEach((r) => {
    const w = (CONTENT.w - totalGap) * (r.w || 1 / regs.length);
    const blocks = r.blocks;
    let est = null;
    for (let adj = 6; adj >= 0; adj--) {
      est = blocks.map((b) => estH(b, w, adj));
      const need = est.reduce((a, e, i) => a + (e ? e.h : MIN_FLEX[blocks[i].type] || 2), 0) + GAP * (blocks.length - 1);
      if (need * 1.06 <= avail || adj === 0) {
        if (need * 1.06 > avail) WARN.push(`#${s.step} region overflow: ${need.toFixed(2)} > ${avail.toFixed(2)}`);
        break;
      }
    }
    const fixed = est.reduce((a, e) => a + (e ? e.h : 0), 0);
    const nFlex = est.filter((e) => !e).length;
    const flexH = nFlex ? (avail - fixed - GAP * (blocks.length - 1)) / nFlex : 0;
    let y = top;
    blocks.forEach((b, i) => {
      if (est[i] && est[i].fs) b._fs = est[i].fs;
      const h = est[i] ? est[i].h : flexH;
      placed.push({ b, x, y, w, h });
      y += h + GAP;
    });
    x += w + GAP;
  });
  s._placed = placed;
  return placed;
}
function scaleCols(colW, w) {
  if (!colW) return null;
  const t = colW.reduce((a, c) => a + c, 0);
  return colW.map((c) => (c * w) / t);
}

// ---------- PPTX renderer ----------
async function buildPptx() {
  const pres = new pptxgen();
  pres.layout = "LAYOUT_WIDE";
  pres.title = "Felodipine ER — full work summary";
  pres.author = "Nikita Ramwani";
  pres.theme = { headFontFace: THEME.headFontFace, bodyFontFace: THEME.bodyFontFace };
  const C = pres.SchemeColor;
  const col = (n) => C[n];

  pres.defineSlideMaster({
    title: "TITLE", background: { color: C.text2 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.0, w: 11.7, h: 1.6, fontFace: THEME.headFontFace, fontSize: 40, bold: true, color: C.background1, valign: "bottom" }, text: "" } },
      { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.8, w: 11.7, h: 1.0, fontSize: 18, color: C.background2, valign: "top" }, text: "" } },
    ],
  });
  pres.defineSlideMaster({
    title: "SECTION", background: { color: C.accent1 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: 0.8, y: 2.4, w: 11.7, h: 1.3, fontFace: THEME.headFontFace, fontSize: 38, bold: true, color: C.background1, valign: "bottom" }, text: "" } },
      { placeholder: { options: { name: "body", type: "body", x: 0.8, y: 3.85, w: 11.7, h: 0.9, fontSize: 18, color: C.background2, valign: "top" }, text: "" } },
    ],
  });
  pres.defineSlideMaster({
    title: "CONTENT", background: { color: C.background1 },
    objects: [
      { placeholder: { options: { name: "title", type: "title", x: 0.6, y: 0.82, w: 12.13, h: 0.75, fontFace: THEME.headFontFace, fontSize: 26, bold: true, color: C.text2, valign: "middle", margin: 0 }, text: "" } },
      { text: { text: "Felodipine ER · MBBE drug 2 · full work summary", options: { x: 0.6, y: 7.05, w: 8, h: 0.3, fontSize: 10, color: C.accent3, margin: 0 } } },
    ],
    slideNumber: { x: 12.2, y: 7.05, w: 0.53, h: 0.3, fontSize: 10, color: C.accent3, align: "right" },
  });

  let secTitle = "Overview";
  pres.addSection({ title: secTitle });
  for (const s of S) {
    if (s.kind === "title") {
      const sl = pres.addSlide({ masterName: "TITLE", sectionTitle: secTitle });
      sl.addText(s.title, { placeholder: "title" });
      sl.addText(s.sub, { placeholder: "body" });
      sl.addText(s.meta, { x: 0.8, y: 5.9, w: 11.7, h: 0.4, fontSize: 14, color: C.background2, isTextBox: true, margin: 0 });
      continue;
    }
    if (s.kind === "section") {
      secTitle = s.title;
      pres.addSection({ title: secTitle });
      const sl = pres.addSlide({ masterName: "SECTION", sectionTitle: secTitle });
      sl.addText(s.title, { placeholder: "title" });
      sl.addText(s.sub, { placeholder: "body" });
      continue;
    }
    const sl = pres.addSlide({ masterName: "CONTENT", sectionTitle: secTitle });
    sl.addText(s.title, { placeholder: "title" });
    checkFit(`#${s.step} title`, [s.title], 12.1, 0.8, 26, 1.15);
    // step badge motif
    sl.addShape(pres.shapes.OVAL, { x: 0.6, y: 0.32, w: 0.42, h: 0.42, fill: { color: C.accent1 }, line: { color: C.accent1 }, objectName: "step badge" });
    sl.addText(String(s.step), { x: 0.6, y: 0.32, w: 0.42, h: 0.42, fontSize: 13, bold: true, color: C.background1, align: "center", valign: "middle", margin: 0, isTextBox: true });
    for (const p of layoutSlide(s)) renderPptxBlock(pres, sl, p, col, s);
    if (s.footnote) {
      sl.addText(s.footnote, { x: 0.6, y: 6.4, w: 12.13, h: 0.5, fontSize: 11, italic: true, color: C.accent3, valign: "top", margin: 0, isTextBox: true });
      checkFit(`#${s.step} footnote`, [s.footnote], 12.1, 0.5, 11);
    }
  }
  const file = path.join(OUT, `${BASENAME}.pptx`);
  await pres.writeFile({ fileName: file });
  return file;
}

function textRuns(item, base) {
  return runs(item).map((r) => ({ text: r.text, options: { ...base, bold: r.bold } }));
}

function renderPptxBlock(pres, sl, p, col, s) {
  const { b, x, y, w, h } = p;
  const label = `#${s.step} ${b.type}`;
  if (b.type === "bullets") {
    const fs = fsOf(b);
    const arr = [];
    b.items.forEach((it, i) => {
      const rr = textRuns(it, { fontSize: fs, color: col("text1") });
      rr[0].options.bullet = { indent: 14 };
      rr[0].options.paraSpaceAfter = 6;
      rr[rr.length - 1].options.breakLine = i < b.items.length - 1;
      arr.push(...rr);
    });
    sl.addText(arr, { x, y, w, h, valign: "top", margin: 0.02, isTextBox: true, fit: "none" });
    checkFit(label, b.items.map(plain), w - 0.3, h, fs, 1.2, 0.09);
  } else if (b.type === "table") {
    const fs = fsOf(b);
    const colW = scaleCols(b.colW, w);
    const hl = typeof b.highlight === "function" ? b.highlight : () => false;
    const head = b.head.map((c) => ({ text: plain(c), options: { bold: true, color: col("background1"), fill: { color: col("text2") }, fontSize: fs } }));
    const body = b.rows.map((r, ri) => r.map((c) => ({
      text: textRuns(c, { fontSize: fs, color: col("text1") }),
      options: { fill: { color: hl(r) ? col("background2") : col("background1") } },
    })));
    sl.addTable([head, ...body], { x, y, w, colW, fontSize: fs, margin: [0.04, 0.06, 0.04, 0.06], border: { type: "solid", pt: 0.5, color: "C9D6D9" }, valign: "middle", autoPage: false });
  } else if (b.type === "stats") {
    const cols = b.cols || b.items.length;
    const rowsN = Math.ceil(b.items.length / cols);
    const cw = (w - GAP * (cols - 1)) / cols, ch = (h - 0.2 * (rowsN - 1)) / rowsN;
    b.items.forEach((it, i) => {
      const cx = x + (i % cols) * (cw + GAP), cy = y + Math.floor(i / cols) * (ch + 0.2);
      sl.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: cx, y: cy, w: cw, h: ch, rectRadius: 0.08, fill: { color: col("background2") }, line: { color: col("background2") }, objectName: "stat card" });
      const vfs = it.v.length > 12 ? 20 : 26;
      sl.addText(it.v, { x: cx + 0.15, y: cy + 0.1, w: cw - 0.3, h: ch * 0.5, fontSize: vfs, bold: true, color: col("accent1"), fontFace: THEME.headFontFace, valign: "bottom", margin: 0, isTextBox: true });
      sl.addText(it.l, { x: cx + 0.15, y: cy + 0.1 + ch * 0.5, w: cw - 0.3, h: ch * 0.5 - 0.15, fontSize: 12, color: col("text1"), valign: "top", margin: 0, isTextBox: true });
      checkFit(label + " value", [it.v], cw - 0.3, ch * 0.5, vfs, 1.1);
      checkFit(label + " label", [it.l], cw - 0.3, ch * 0.5 - 0.15, 12);
    });
  } else if (b.type === "callout") {
    sl.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08, fill: { color: col("text2") }, line: { color: col("text2") }, objectName: "callout" });
    sl.addText(textRuns(b.text, { fontSize: 15, color: col("background1") }), { x: x + 0.2, y, w: w - 0.4, h, valign: "middle", margin: 0, isTextBox: true });
    checkFit(label, [b.text], w - 0.4, h - 0.1, 15);
  } else if (b.type === "note") {
    sl.addText(b.text, { x, y, w, h, fontSize: 11, italic: true, color: col("accent3"), valign: "top", margin: 0, isTextBox: true });
  } else if (b.type === "cards") {
    const tone = b.tone || "accent1";
    sl.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y, w, h, rectRadius: 0.08, fill: { color: col("background2") }, line: { color: col("background2") }, objectName: "card" });
    sl.addShape(pres.shapes.OVAL, { x: x + 0.2, y: y + 0.22, w: 0.22, h: 0.22, fill: { color: col(tone) }, line: { color: col(tone) }, objectName: "card dot" });
    sl.addText(b.title, { x: x + 0.55, y: y + 0.12, w: w - 0.75, h: 0.42, fontSize: 15, bold: true, color: col("text2"), valign: "middle", margin: 0, isTextBox: true });
    const fs = fsOf(b);
    const arr = [];
    b.items.forEach((it, i) => {
      const rr = textRuns(it, { fontSize: fs, color: col("text1") });
      rr[0].options.bullet = { indent: 14 };
      rr[0].options.paraSpaceAfter = 6;
      rr[rr.length - 1].options.breakLine = i < b.items.length - 1;
      arr.push(...rr);
    });
    sl.addText(arr, { x: x + 0.2, y: y + 0.65, w: w - 0.4, h: h - 0.8, valign: "top", margin: 0, isTextBox: true });
    checkFit(label, b.items.map(plain), w - 0.7, h - 0.8, fs, 1.2, 0.09);
  } else if (b.type === "image") {
    const f = img(b.src), d = pngSize(f);
    const r = Math.min(w / d.w, h / d.h);
    const iw = d.w * r, ih = d.h * r;
    sl.addImage({ path: f, x: x + (w - iw) / 2, y: y + (h - ih) / 2, w: iw, h: ih, altText: path.basename(b.src) });
  } else if (b.type === "flow") {
    const n = b.steps.length;
    if (b.vertical) {
      const sh = (h - 0.2 * (n - 1)) / n;
      b.steps.forEach((st, i) => {
        const cy = y + i * (sh + 0.2);
        sl.addShape(pres.shapes.ROUNDED_RECTANGLE, { x, y: cy, w, h: sh, rectRadius: 0.08, fill: { color: col("background2") }, line: { color: col("background2") }, objectName: "flow step" });
        sl.addShape(pres.shapes.OVAL, { x: x + 0.18, y: cy + (sh - 0.5) / 2, w: 0.5, h: 0.5, fill: { color: col("accent1") }, line: { color: col("accent1") }, objectName: "flow badge" });
        sl.addText(String(i + 1), { x: x + 0.18, y: cy + (sh - 0.5) / 2, w: 0.5, h: 0.5, fontSize: 14, bold: true, color: col("background1"), align: "center", valign: "middle", margin: 0, isTextBox: true });
        sl.addText([{ text: st.h, options: { bold: true, fontSize: 16, color: col("text2"), breakLine: true } }, { text: st.t, options: { fontSize: 15, color: col("text1") } }],
          { x: x + 0.85, y: cy, w: w - 1.0, h: sh, valign: "middle", margin: 0, isTextBox: true });
        checkFit(label, [st.h, st.t], w - 1.0, sh, 13.5);
      });
    } else {
      const aw = 0.25, bw = (w - aw * (n - 1)) / n;
      b.steps.forEach((st, i) => {
        const cx = x + i * (bw + aw);
        sl.addShape(pres.shapes.ROUNDED_RECTANGLE, { x: cx, y, w: bw, h, rectRadius: 0.08, fill: { color: i === n - 1 ? col("accent2") : col("text2") }, line: { color: i === n - 1 ? col("accent2") : col("text2") }, objectName: "flow step" });
        sl.addText([{ text: st.h, options: { bold: true, fontSize: 16, color: col("background1"), breakLine: true } }, { text: st.t, options: { fontSize: 14, color: col("background1") } }],
          { x: cx + 0.15, y: y + 0.1, w: bw - 0.3, h: h - 0.2, valign: "top", margin: 0, isTextBox: true });
        checkFit(label, [st.h, st.t], bw - 0.3, h - 0.2, 14);
        if (i < n - 1) sl.addShape(pres.shapes.RIGHT_ARROW, { x: cx + bw + 0.04, y: y + h / 2 - 0.1, w: 0.17, h: 0.2, fill: { color: col("accent3") }, line: { color: col("accent3") }, objectName: "flow arrow" });
      });
    }
  } else if (b.type === "chart") {
    const hex = [THEME.colors.accent1, THEME.colors.accent2, THEME.colors.accent3, THEME.colors.accent4];
    const data = b.series.map((se) => ({ name: se.name, labels: b.labels, values: se.values }));
    const isLine = b.kind === "line";
    const opts = {
      x, y, w, h,
      showTitle: true, title: b.title, titleFontSize: 14, titleColor: THEME.colors.dk2, titleFontFace: "+mn-lt",
      catAxisLabelColor: THEME.colors.dk1, valAxisLabelColor: THEME.colors.dk1, catAxisLabelFontSize: 11, valAxisLabelFontSize: 11,
      catAxisLabelFontFace: "+mn-lt", valAxisLabelFontFace: "+mn-lt",
      valGridLine: { color: "D5DEE0", size: 0.5 }, catGridLine: { style: "none" },
      valAxisMinVal: b.yMin, valAxisMaxVal: b.yMax, valAxisLabelFormatCode: "0.0",
      showLegend: !!b.legend || isLine, legendPos: "b", legendFontSize: 11, legendFontFace: "+mn-lt",
      showValue: !isLine, dataLabelPosition: "outEnd", dataLabelFontSize: 11, dataLabelFormatCode: b.valueFmt || "0.00", dataLabelColor: THEME.colors.dk1, dataLabelFontFace: "+mn-lt",
      chartColors: isLine ? [THEME.colors.accent1, THEME.colors.accent4, THEME.colors.accent4] : hex.slice(0, b.series.length),
    };
    if (b.xTitle) Object.assign(opts, { showCatAxisTitle: true, catAxisTitle: b.xTitle, catAxisTitleFontSize: 11, catAxisTitleColor: THEME.colors.dk1, catAxisTitleFontFace: "+mn-lt" });
    if (isLine) {
      // pptxgenjs writes one dash style for every series; fixLineCharts() sets the limit lines dashed afterwards.
      Object.assign(opts, { lineSize: 2, lineDataSymbol: "circle", lineDataSymbolSize: 7, lineDash: "solid" });
      sl.addChart(pres.charts.LINE, data, opts);
    } else {
      Object.assign(opts, { barDir: "col", barGapWidthPct: 60 });
      sl.addChart(pres.charts.BAR, data, opts);
    }
  }
}

// ---------- HTML renderer ----------
const PX = 96; // px per inch on a 1280×720 slide
const pt = (v) => (v * 96) / 72;
function htmlRuns(s) {
  return runs(s).map((r) => (r.bold ? `<b>${esc(r.text)}</b>` : esc(r.text))).join("");
}
function box(x, y, w, h, inner, cls = "", style = "") {
  return `<div class="${cls}" data-h="${Math.round(h * PX)}" style="left:${x * PX}px;top:${y * PX}px;width:${w * PX}px;height:${h * PX}px;${style}">${inner}</div>`;
}
function svgChart(b, w, h) {
  const Wp = w * PX, Hp = h * PX;
  const padL = 46, padR = 14, padT = 34, padB = (b.xTitle ? 58 : 40) + (b.legend || b.kind === "line" ? 22 : 0);
  const iw = Wp - padL - padR, ih = Hp - padT - padB;
  const yv = (v) => padT + ih - ((v - b.yMin) / (b.yMax - b.yMin)) * ih;
  const n = b.labels.length, step = iw / n;
  const colors = b.kind === "line" ? [HEX.accent1, HEX.accent4, HEX.accent4] : [HEX.accent1, HEX.accent2, HEX.accent3];
  let g = `<text x="${Wp / 2}" y="20" class="ct">${esc(b.title)}</text>`;
  const ticks = 5;
  for (let i = 0; i <= ticks; i++) {
    const v = b.yMin + ((b.yMax - b.yMin) * i) / ticks, yy = yv(v);
    g += `<line x1="${padL}" x2="${padL + iw}" y1="${yy}" y2="${yy}" class="grid"/><text x="${padL - 6}" y="${yy + 4}" class="ax" text-anchor="end">${v.toFixed(1)}</text>`;
  }
  b.labels.forEach((l, i) => { g += `<text x="${padL + step * (i + 0.5)}" y="${padT + ih + 18}" class="ax" text-anchor="middle">${esc(l)}</text>`; });
  if (b.xTitle) g += `<text x="${padL + iw / 2}" y="${padT + ih + 40}" class="ax" text-anchor="middle">${esc(b.xTitle)}</text>`;
  if (b.kind === "line") {
    b.series.forEach((se, si) => {
      const pts = se.values.map((v, i) => `${padL + step * (i + 0.5)},${yv(v)}`).join(" ");
      g += `<polyline points="${pts}" fill="none" stroke="#${colors[si]}" stroke-width="${se.dash ? 1.5 : 2.5}" ${se.dash ? 'stroke-dasharray="6 4"' : ""}/>`;
      if (!se.dash) se.values.forEach((v, i) => { g += `<circle cx="${padL + step * (i + 0.5)}" cy="${yv(v)}" r="4" fill="#${colors[si]}"><title>${esc(b.labels[i])}: ${v}</title></circle>`; });
    });
  } else {
    const k = b.series.length, bw = (step * 0.62) / k;
    b.series.forEach((se, si) => se.values.forEach((v, i) => {
      const bx = padL + step * i + step * 0.19 + si * bw, by = yv(v);
      g += `<rect x="${bx}" y="${by}" width="${bw - 3}" height="${padT + ih - by}" rx="2" fill="#${colors[si]}"><title>${esc(se.name)} · ${esc(b.labels[i])}: ${v}</title></rect>`;
      g += `<text x="${bx + (bw - 3) / 2}" y="${by - 5}" class="ax" text-anchor="middle">${v.toFixed(b.valueFmt === "0.0" ? 1 : 2)}</text>`;
    }));
  }
  if (b.legend || b.kind === "line") {
    let lx = padL;
    const ly = Hp - 10;
    b.series.forEach((se, si) => {
      g += `<line x1="${lx}" x2="${lx + 18}" y1="${ly - 4}" y2="${ly - 4}" stroke="#${colors[si]}" stroke-width="3" ${se.dash ? 'stroke-dasharray="5 3"' : ""}/><text x="${lx + 24}" y="${ly}" class="ax">${esc(se.name)}</text>`;
      lx += 40 + se.name.length * 7;
    });
  }
  return `<svg viewBox="0 0 ${Wp} ${Hp}" width="${Wp}" height="${Hp}" role="img" aria-label="${esc(b.title)}">${g}</svg>`;
}
function renderHtmlBlock(p) {
  const { b, x, y, w, h } = p;
  if (b.type === "bullets") {
    const fs = fsOf(b);
    return box(x, y, w, h, `<ul style="font-size:${pt(fs)}px">${b.items.map((i) => `<li>${htmlRuns(i)}</li>`).join("")}</ul>`, "blk");
  }
  if (b.type === "table") {
    const fs = fsOf(b), colW = scaleCols(b.colW, w);
    const hl = typeof b.highlight === "function" ? b.highlight : () => false;
    const cg = colW ? `<colgroup>${colW.map((c) => `<col style="width:${c * PX}px">`).join("")}</colgroup>` : "";
    const head = `<tr>${b.head.map((c) => `<th>${esc(plain(c))}</th>`).join("")}</tr>`;
    const body = b.rows.map((r) => `<tr class="${hl(r) ? "hl" : ""}">${r.map((c) => `<td>${htmlRuns(c)}</td>`).join("")}</tr>`).join("");
    return box(x, y, w, h, `<table style="font-size:${pt(fs)}px">${cg}${head}${body}</table>`, "blk", "height:auto");
  }
  if (b.type === "stats") {
    const cols = b.cols || b.items.length;
    const rowsN = Math.ceil(b.items.length / cols);
    const cw = (w - GAP * (cols - 1)) / cols, ch = (h - 0.2 * (rowsN - 1)) / rowsN;
    return b.items.map((it, i) => {
      const cx = x + (i % cols) * (cw + GAP), cy = y + Math.floor(i / cols) * (ch + 0.2);
      const vfs = it.v.length > 12 ? 20 : 26;
      return box(cx, cy, cw, ch, `<div class="sv" style="font-size:${pt(vfs)}px;height:50%">${esc(it.v)}</div><div class="sl">${esc(it.l)}</div>`, "stat");
    }).join("");
  }
  if (b.type === "callout") return box(x, y, w, h, `<div>${htmlRuns(b.text)}</div>`, "callout");
  if (b.type === "note") return box(x, y, w, h, esc(b.text), "note");
  if (b.type === "cards") {
    return box(x, y, w, h, `<div class="ct2"><span class="dot" style="background:var(--${b.tone || "accent1"})"></span>${esc(b.title)}</div><ul style="font-size:${pt(fsOf(b))}px">${b.items.map((i) => `<li>${htmlRuns(i)}</li>`).join("")}</ul>`, "card");
  }
  if (b.type === "image") {
    const f = img(b.src), d = pngSize(f);
    const r = Math.min(w / d.w, h / d.h), iw = d.w * r, ih = d.h * r;
    const data = fs.readFileSync(f).toString("base64");
    return box(x + (w - iw) / 2, y + (h - ih) / 2, iw, ih, `<img alt="${esc(path.basename(b.src))}" src="data:image/png;base64,${data}">`, "blk imgbox");
  }
  if (b.type === "flow") {
    const n = b.steps.length;
    if (b.vertical) {
      const sh = (h - 0.2 * (n - 1)) / n;
      return b.steps.map((st, i) => box(x, y + i * (sh + 0.2), w, sh, `<span class="badge">${i + 1}</span><div><b>${esc(st.h)}</b><br>${esc(st.t)}</div>`, "vflow")).join("");
    }
    const aw = 0.25, bw = (w - aw * (n - 1)) / n;
    return b.steps.map((st, i) => box(x + i * (bw + aw), y, bw, h, `<b>${esc(st.h)}</b><br>${esc(st.t)}`, "hflow" + (i === n - 1 ? " last" : "")) +
      (i < n - 1 ? box(x + i * (bw + aw) + bw, y + h / 2 - 0.12, aw, 0.24, "▶", "arrow") : "")).join("");
  }
  if (b.type === "chart") return box(x, y, w, h, svgChart(b, w, h), "blk");
  return "";
}
function buildHtml() {
  const slidesHtml = S.map((s, i) => {
    if (s.kind === "title") return `<section class="slide dark" data-i="${i}">${box(0.8, 2.0, 11.7, 1.6, `<h1>${esc(s.title)}</h1>`, "blk tb")}${box(0.8, 3.8, 11.7, 1.0, esc(s.sub), "blk sub")}${box(0.8, 5.9, 11.7, 0.4, esc(s.meta), "blk meta")}</section>`;
    if (s.kind === "section") return `<section class="slide sec" data-i="${i}">${box(0.8, 2.4, 11.7, 1.3, `<h1>${esc(s.title)}</h1>`, "blk tb")}${box(0.8, 3.85, 11.7, 0.9, esc(s.sub), "blk sub")}</section>`;
    const parts = layoutSlide(s).map(renderHtmlBlock).join("");
    return `<section class="slide" data-i="${i}">` +
      box(0.6, 0.32, 0.42, 0.42, String(s.step), "badge0") +
      box(0.6, 0.82, 12.13, 0.75, `<h2>${esc(s.title)}</h2>`, "blk tt") + parts +
      (s.footnote ? box(0.6, 6.4, 12.13, 0.5, esc(s.footnote), "note") : "") +
      box(0.6, 7.05, 8, 0.3, "Felodipine ER · MBBE drug 2 · full work summary", "foot") +
      box(12.2, 7.05, 0.53, 0.3, String(i + 1), "foot num") + `</section>`;
  }).join("\n");
  const c = THEME.colors;
  const html = `<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Felodipine ER Work Summary</title>
<style>
:root{--text1:#${c.dk1};--text2:#${c.dk2};--bg1:#${c.lt1};--bg2:#${c.lt2};--accent1:#${c.accent1};--accent2:#${c.accent2};--accent3:#${c.accent3};--accent4:#${c.accent4};--accent5:#${c.accent5};--page:#d9e2e4}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--page:#11181c}}
:root[data-theme="dark"]{--page:#11181c}
*{box-sizing:border-box}
body{margin:0;background:var(--page);font-family:Calibri,Carlito,"Segoe UI",Arial,sans-serif;color:var(--text1)}
.bar{position:sticky;top:0;z-index:5;display:flex;gap:12px;align-items:center;justify-content:center;padding:10px 16px;background:var(--text2);color:#fff;font-size:14px;flex-wrap:wrap}
.bar button{background:var(--accent1);color:#fff;border:0;border-radius:6px;padding:6px 12px;font:inherit;cursor:pointer}
.deck{display:flex;flex-direction:column;align-items:center;gap:20px;padding:20px 16px 60px}
.frame{width:100%;max-width:1280px;aspect-ratio:16/9;position:relative;overflow:hidden;box-shadow:0 2px 12px rgba(0,0,0,.18);border-radius:4px;background:#fff}
.slide{position:absolute;left:0;top:0;width:1280px;height:720px;transform-origin:0 0;background:var(--bg1);overflow:hidden}
.slide>div{position:absolute}
.dark{background:var(--text2);color:#fff}.sec{background:var(--accent1);color:#fff}
h1{font-family:Cambria,Caladea,Georgia,serif;font-size:${pt(40)}px;margin:0;line-height:1.1;display:flex;align-items:flex-end;height:100%}
.sec h1{font-size:${pt(38)}px}
.tb{display:flex}.sub{font-size:${pt(18)}px;color:var(--bg2);line-height:1.3}.meta{font-size:${pt(14)}px;color:var(--bg2)}
h2{font-family:Cambria,Caladea,Georgia,serif;font-size:${pt(26)}px;color:var(--text2);margin:0;height:100%;display:flex;align-items:center;line-height:1.1}
.badge0{background:var(--accent1);color:#fff;border-radius:50%;font-weight:700;font-size:${pt(13)}px;display:flex;align-items:center;justify-content:center}
.kicker{color:var(--accent2);font-weight:700;font-size:${pt(12)}px;letter-spacing:2px;display:flex;align-items:center}
ul{margin:0;padding-left:20px;line-height:1.2}li{margin-bottom:7px}
table{border-collapse:collapse;width:100%;table-layout:fixed;line-height:1.2}
th{background:var(--text2);color:#fff;text-align:left;font-weight:700}
th,td{border:0.5px solid #C9D6D9;padding:4px 6px;vertical-align:middle}
tr.hl td{background:var(--bg2)}
.stat{background:var(--bg2);border-radius:8px;padding:10px 14px}
.sv{font-family:Cambria,Caladea,Georgia,serif;font-weight:700;color:var(--accent1);display:flex;align-items:flex-end;line-height:1.05}
.sl{font-size:${pt(12)}px;margin-top:4px;line-height:1.2}
.callout{background:var(--text2);color:#fff;border-radius:8px;padding:0 19px;display:flex;align-items:center;font-size:${pt(15)}px;line-height:1.25}
.note{font-size:${pt(11)}px;font-style:italic;color:var(--accent3);line-height:1.25}
.card{background:var(--bg2);border-radius:8px;padding:12px 19px}
.card .ct2{font-weight:700;font-size:${pt(15)}px;color:var(--text2);display:flex;align-items:center;gap:12px;margin-bottom:14px}
.card ul{font-size:${pt(14)}px}
.dot{width:21px;height:21px;border-radius:50%;flex:none}
.imgbox img{width:100%;height:100%;display:block}
.blk svg{display:block}
.vflow{background:var(--bg2);border-radius:8px;display:flex;align-items:center;gap:16px;padding:0 14px 0 17px;font-size:${pt(15)}px;line-height:1.25}
.vflow b{font-size:${pt(16)}px;color:var(--text2)}
.vflow .badge{width:48px;height:48px;border-radius:50%;background:var(--accent1);color:#fff;font-weight:700;display:flex;align-items:center;justify-content:center;flex:none;font-size:${pt(14)}px}
.hflow{background:var(--text2);color:#fff;border-radius:8px;padding:10px 14px;font-size:${pt(14)}px;line-height:1.25}
.hflow b{font-size:${pt(16)}px}.hflow.last{background:var(--accent2)}
.arrow{color:var(--accent3);display:flex;align-items:center;justify-content:center;font-size:14px}
.foot{font-size:${pt(10)}px;color:var(--accent3)}.num{text-align:right}
svg .ct{font:700 ${pt(14)}px Calibri,Carlito,Arial,sans-serif;fill:var(--text2);text-anchor:middle}
svg .ax{font:${pt(11)}px Calibri,Carlito,Arial,sans-serif;fill:var(--text1)}
svg .grid{stroke:#D5DEE0;stroke-width:1}
body.present{overflow:hidden;background:#000}
body.present .bar{display:none}
body.present .deck{padding:0;height:100vh;justify-content:center}
body.present .frame{display:none;max-width:none;width:min(100vw,177.78vh);border-radius:0;box-shadow:none}
body.present .frame.on{display:block}
@page{size:13.333in 7.5in;margin:0}
@media print{.bar{display:none}.deck{padding:0;gap:0}.frame{width:1280px;max-width:1280px;height:720px;box-shadow:none;border-radius:0;page-break-after:always;break-after:page}.slide{transform:none!important}body{background:#fff}}
</style></head><body>
<div class="bar"><span>Felodipine ER — full work summary · ${S.length} slides</span><button id="present">Present (arrow keys, Esc to exit)</button><button onclick="window.print()">Print / save as PDF</button></div>
<main class="deck">
${slidesHtml.split("\n").map((sh) => `<div class="frame">${sh}</div>`).join("\n")}
</main>
<script>
const frames=[...document.querySelectorAll('.frame')];let cur=0;
function fit(){frames.forEach(f=>{const s=f.firstElementChild;if(s)s.style.transform='scale('+(f.clientWidth/1280)+')';});}
addEventListener('resize',fit);fit();
function show(i){cur=Math.max(0,Math.min(frames.length-1,i));frames.forEach((f,k)=>f.classList.toggle('on',k===cur));fit();}
document.getElementById('present').onclick=()=>{document.body.classList.add('present');show(cur);try{document.documentElement.requestFullscreen()}catch(e){}};
addEventListener('keydown',e=>{if(!document.body.classList.contains('present'))return;
 if(['ArrowRight','PageDown',' '].includes(e.key)){show(cur+1);e.preventDefault();}
 else if(['ArrowLeft','PageUp'].includes(e.key)){show(cur-1);e.preventDefault();}
 else if(e.key==='Escape'){document.body.classList.remove('present');frames.forEach(f=>f.classList.remove('on'));fit();}});
document.addEventListener('fullscreenchange',()=>{if(!document.fullscreenElement){document.body.classList.remove('present');frames.forEach(f=>f.classList.remove('on'));fit();}});
</script></body></html>`;
  const file = path.join(OUT, `${BASENAME}.html`);
  fs.writeFileSync(file, html);
  return file;
}


// pptxgenjs line-chart XML needs three fixes for PowerPoint: dashed limit series (one dash
// style is written for all series), no dLblPos (outEnd is not valid on a line chart), and
// exactly the two axis ids the chart declares.
async function fixLineCharts(file) {
  const JSZip = require("jszip");
  const zip = await JSZip.loadAsync(fs.readFileSync(file));
  for (const name of Object.keys(zip.files).filter((n) => /^ppt\/charts\/chart\d+\.xml$/.test(n))) {
    let x = await zip.file(name).async("string");
    if (!x.includes("<c:lineChart>")) continue;
    const declared = new Set([...x.matchAll(/<c:(?:catAx|valAx)>\s*<c:axId val="(\d+)"/g)].map((m) => m[1]));
    x = x.replace(/<c:lineChart>[\s\S]*?<\/c:lineChart>/g, (g) => g
      .replace(/<c:dLblPos val="[^"]*"\/>/g, "")
      .replace(/<c:axId val="(\d+)"\/>/g, (m, id) => (declared.has(id) ? m : ""))
      .replace(/<c:ser>[\s\S]*?<\/c:ser>/g, (ser) => {
        const nm = (ser.match(/<c:v>([^<]*)<\/c:v>/) || [])[1];
        if (!LIMIT_SERIES.has(nm)) return ser;
        return ser.replace(/<a:prstDash val="[^"]*"\/>/g, '<a:prstDash val="dash"/>')
          .replace(/<a:ln w="\d+"/, '<a:ln w="15875"')
          .replace(/<c:marker>[\s\S]*?<\/c:marker>/, '<c:marker><c:symbol val="none"/></c:marker>');
      }));
    zip.file(name, x);
  }
  fs.writeFileSync(file, await zip.generateAsync({ type: "nodebuffer", compression: "DEFLATE" }));
}
const LIMIT_SERIES = new Set(["1.25", "0.80"]);

(async () => {
  const p = await buildPptx();
  const { applyTheme } = require(process.env.APPLY_THEME);
  await applyTheme(p, THEME);
  await fixLineCharts(p);
  const h = buildHtml();
  console.log("slides:", S.length);
  console.log(p);
  console.log(h);
  if (WARN.length) console.log("FIT WARNINGS:\n" + WARN.join("\n"));
})();
