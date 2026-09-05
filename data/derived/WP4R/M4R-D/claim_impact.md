# M4R-D claim impact (ROUND1, 2026-09-05)

## CLM08 (composition-coupled turnover partitions architecture, differs from adaptive sites)
- D1 branch-aware: SUPPORTS + STRENGTHENS. Coupled-site branch substitutions follow branch
  GC/AT direction in 250/316 informative events (0.791, CI 0.743–0.832, binom p=2.4e-26);
  uncoupled 5538/8821 (0.628); Fisher coupled-vs-uncoupled OR=2.25, p=9.5e-10.
  Multi-nt substitutions more concordant (0.667) than 1-nt (0.546). Terminal correlation
  retained as descriptive layer; branch layer is now the stronger evidence.
- D2 architecture: MIXED but informative. Project-convention LCR (window64/ent1.5) flags 0%
  of core-OG sites in both groups — convention uninformative on proteins (disclosed limitation,
  not biology). Homopolymer-run/polyN enrichment of coupled sites (2/223 vs 5/9662, OR=17.5,
  p=0.0099; counts tiny). Honest headline: **within conserved single-copy core genes, ~99% of
  composition-coupled turnover falls in NON-LCR sequence** — i.e. churn is not repeat-confined
  WITHIN the core-gene sampling frame. **SAMPLING-FRAME CAVEAT (adversarial review 2026-09-05):
  core-181 construction excludes LCR-heavy/fast-evolving proteins (D6 found 0 AP2/PUF among
  core-181), so "churn is not in LCRs" is partly built into the frame; whole-proteome
  architecture mapping (D2b, unbounded by core OGs) is a mandatory round-2 item before any
  Figure-3-level architecture claim.**
- D3 conservation-matched controls (conservation-only matching; domain/essentiality/structure
  matching NOT yet done — round-2 item): SUPPORTS weak version. 14/14 known resistance sites uncoupled
  (r undefined); conservation-matched backgrounds: coupling 0% at cons=1.0 sites vs 9% at
  mid-cons sites; adaptive-sub Grantham spans full coupled distribution (quantiles 0.0–1.0,
  median 74.5 vs 77.3). Keep: "canonical resistance mutations occur at highly constrained
  sites, whereas composition-driven turnover dominates a different part of sequence space."
  NO dual-adaptive-axes claim.
- D6 natural experiment: SUGGESTIVE pilot only (4 mafft-pilot AP2/PUF OGs + 5 CHROM core OGs;
  pole Asn diffs +0.0–+2.4pp, anchors 53–93% conserved). AP2/PUF absent from core-181 by
  construction (multi-copy/fast-evolving) — generality limit disclosed.

## CLM09 (regulatory Asn/LCR/IDR enrichment)
- SUPPORTS for AP2/PUF; REJECTS for CHROM. AP2: crude OR=34.2 (p=1.6e-17);
  length×repeat-stratified MH OR=27.4, CMH p=5.4e-11, positive in 2/2 informative strata
  (LxTract5+: regulator median Asn 0.249 vs bg 0.163, MW p=6.3e-12). PUF/RNA:
  MH OR=3.03, p=2.6e-3. CHROM: null (OR=1.19, p=0.83).
  Counterfactual holds: given comparable length and repeat content, AP2/PUF proteins are
  unusually Asn-rich — not reducible to "more IDR".
- Constraint context (Zhang Pf S5): AP2s tolerate insertions (no-insertion frac 0.39 vs 0.62;
  MIS 1.00 vs 0.39, p=1.7e-4) — not blood-stage essential; lifecycle-conditional
  essentiality is the testable frame, not generic essentiality.

## CLM03/CLM01 (supporting)
- Branch-aware concordance is independent support that composition propagates into AA space
  along phylogeny (CLM03), built on established multi-state transitions (CLM01).

## Gate impact
- D-backbone verdict: CONDITIONAL (strong components D1/D3/D5; D2 honest-mixed; D6 pilot).
  Figure 2–3 candidacy holds. Round 2 needs: formal D6 (≥2 independent transitions),
  Pf8 low-frequency boundary, official-OG sensitivity.

## Limitations
- kmerRBH core-181 ascertainment (AP2/PUF absent; EXT drug-gene OGs outside v2 stats).
- Fitch parsimony (646 ambiguous-root sites flagged, retained — exclusion changes nothing
  directionally; uncertainty disclosed in aux table).
- Convention LCR caller uninformative on proteins (reported, sensitivity thr2.0 shown).
- Association only — no mechanism claimed.

## ROUND2 append (M4R-D-ROUND2)
- D2b verdict: SURVIVES (frame-internal). Coupled-harboring Pf genes (91) are NOT more LCR-rich
  than length/age-matched background (lcr med 0.0 vs 0.0, MW p=0.10; strat-diff 0.0); they are
  Asn-DEPLETED (0.064 vs 0.119, p=3.9e-33) with shorter polyN (2 vs 3, p=2.4e-18) and MORE
  InterPro domains (5 vs 1, p=5.7e-28) — structured conserved proteins, not repeat proteins.
  Frame bias quantified: core-181 members shorter (259 vs 472 aa, p=7.8e-23) and Asn-poorer
  (0.056 vs 0.120, p=4.9e-75) than non-core; "churn not repeat-confined" holds WITHIN frame.
- D6 formal: 2/3 transitions composition-consistent. T2_GC_vivax: AP2 4/4 (p=0.0625),
  PUF_RNA 21/27 (p=0.0030); T3_rodent: AP2 3/4, PUF_RNA 18/25 (p=0.022); T1_AT_Lav negative
  (AP2 2/4, PUF 9/25; registered polarity tiny: 0.19 vs 0.18 GC — weak transition by construction).
  CHROM null flat on all 3. Support rule (>=2 independent) MET by T2+T3.
- Branch-LOO: coupled-vs-uncoupled Fisher OR range 1.71 (vivax-drop, p=0.0065) to 4.35
  (piroplasm-drop); significant in all 7. Weakest under vivax-drop (high-GC pole removed).
  Piroplasm-drop threshold jumps to 0.142 (dGC compression) — disclosed.
- D3 multidim (partial upgrade): domain-shared bg where available — K13/Kelch 531 sites @4.5%
  coupled vs known 0%; MDR1/ABC 5056 sites @2.0% vs known 0% (broad domain, weak restriction);
  CRT/DHFR/DHPS have no domain-shared core sites (family absent from core set).
  Within-OG control structurally impossible (all 14 known-gene OGs absent from core alignments —
  fast-evolving drug genes fail QC; consistent with D6 frame finding). Zhang MIS parsed
  (K13/CRT/DHFR/MDR1 all 0.12–0.13, dispensable-bin); matched-ess analysis deferred as
  uninformative (no MIS contrast among the 5 genes). Scope renamed: conservation-matched +
  domain-shared where available.
- Pf8 boundary: 34/34 markers segregating (call rate ~1.0); crt76 derived 0.56, dhfr51 0.83,
  dhfr108 0.95 (strict single-clone subset consistent). Adaptation segregates in populations
  at cross-species-conserved sites — D3 separation holds at population layer.
  CNV: GCH1 0.28, PM2/PM3 0.28, MDR1 0.20, CRT 0.086 amplifications; HRP2/3 deletions ~0.45.
  SNP-level burden DEFERRED (Zarr streaming beyond round-2 budget).
- Prohibited extrapolations respected: no mechanism language; T1 negative kept (not rescued);
  D2b claim kept frame-internal; no new downloads.
