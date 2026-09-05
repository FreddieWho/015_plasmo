# M4R-D ROUND1 REPORT (2026-09-05)

**Verdict: CONDITIONAL** (strong D1/D3/D5 components; D2 honest-mixed; D6 pilot-only).
D can carry the evolutionary backbone (Figure 2–3 candidacy holds) pending round-2
formalization. Association only — no mechanism claimed.

## Inputs (read-only, no downloads)
- L2v2 site table (29,558 parsimony sites; 223 coupled q<0.05 / 9,662 uncoupled q>=0.1),
  adaptive-sub Grantham (14), alignment QC (175 pass), ML supermatrix tree (16 spp),
  181 core + 5 EXT alignments; M1 GC states; Grantham refs; NCBI + PlasmoDB-71 Pf proteome;
  Elsworth S3a (5,266 Pk), Zhang S5 (5,403 Pf), Oberstaller supplement (in hand, unused);
  L1d keyword lists. Scripts: run_m4rd.py / round2 / round3 (seed 20260905; sc env).

## D1 branch-aware L2 — STRONG, upgrades terminal correlation
- 24,520 branch substitutions mapped (Fitch + squared-change GC; 646 ambiguous-root flagged).
- Informative branches (|dGC|>=median 0.037): overall concordance 7,927/13,604 = **0.583**
  (CI 0.574–0.591, p=3.2e-83).
- **Coupled sites: 250/316 = 0.791** (CI 0.743–0.832, p=2.4e-26) vs uncoupled 0.628;
  Fisher OR=**2.25**, p=9.5e-10.
- Multi-nt subs 0.667 vs 1-nt 0.546 (accessibility gradient consistent with composition-driven change).
- Effect: branch-aware layer now stronger than terminal correlation (which is retained descriptive).

## D2 architecture partition — MIXED, honest headline (sampling-frame caveat added post-review)
- Convention LCR (window64/ent1.5): 0.0% both groups — caller uninformative on proteins
  (convention limitation, disclosed; thr2.5 saturates at 99%, thr2.0: 1/223 vs 2/9662, OR=21.8, p=0.066).
- Homopolymer≥5: 2/223 vs 5/9662, OR=**17.5**, p=0.0099 (tiny counts — weak, not load-bearing).
- Headline: **within the conserved single-copy core-gene frame, ~99% of coupled turnover is
  non-LCR** with radical chemistry (median Grantham 77.3 vs 66.9). Caveat: core-181 construction
  excludes LCR-heavy/fast-evolving proteins (0 AP2/PUF in core-181), so "not repeat-confined"
  is partly frame-built. Whole-proteome architecture mapping (D2b) required before Figure-3 claims.
  Logistic: conservation coef −1.31 (coupling avoids the most conserved columns), n_class ≈ 0.

## D3 conservation-matched adaptive-site controls — STRONG descriptive, weak claim kept
- 14/14 known resistance positions uncoupled (r undefined; cons 0.50–1.00).
- Cons-matched backgrounds: coupling 0% at cons=1.0 (n=6,454) vs 9.1% at cons=0.5;
  bg median |r| 0.27→0.40, Grantham 45→74.6, n_aa 2→4 along falling conservation.
- Adaptive subs span the full coupled Grantham range (quantiles 0.00–1.00; median 74.5 vs 77.3).
- Verdict: separation is real but conservation-graded → keep weak wording only
  (constrained resistance sites vs turnover elsewhere). No dual-axes claim.

## D5 regulatory enrichment + counterfactual — STRONG (AP2/PUF), null (CHROM)
- AP2 (n=28): crude OR=34.2, p=1.6e-17; length×repeat-stratified **MH OR=27.4, CMH p=5.4e-11**,
  positive in 2/2 strata (large+tract stratum: 0.249 vs 0.163, MW p=6.3e-12).
- PUF/RNA (n=56): MH OR=3.03, p=2.6e-3. CHROM (n=56): null (OR=1.19, p=0.83).
- Counterfactual holds: given comparable length+repeat content, AP2/PUF IDRs/proteins are
  unusually Asn-rich.
- Constraint: AP2s tolerate blood-stage insertions (Zhang: no-insertion 0.39 vs 0.62;
  MIS 1.00 vs 0.39, p=1.7e-4) → lifecycle-conditional frame, not generic essentiality.
  Pk background: 2,037 essential / 2,124 dispensable / 1,105 intermediate.

## D6 natural experiment — SUGGESTIVE pilot
- Structural limit found: **0 AP2/PUF genes in core-181** (multi-copy/fast-evolving, excluded
  by construction) → D6 cannot run on L2 alignments for the key families.
- Salvage (mafft pilot, 2 AP2 + 2 PUF wide-coverage OGs + 5 CHROM core OGs): domain anchors
  53–93% conserved; low-GC pole Asn exceeds high-GC pole by +0.9–+2.4pp in 4/4 pilot OGs,
  ~0 in CHROM. Directionally consistent, underpowered — UNRESOLVED, needs round-2 formal test
  with ≥2 independent transitions.

## D4 mutation→fixation — boundary only
- No resistant×starvation proteome exists (GAP-1 stands). Delivered: Zhang/Elsworth constraint
  context + Hamilton MA cited as single-species control. Pf8 low-frequency boundary deferred
  to round 2.

## What would change the verdict
- →STRONG: round-2 D6 formalized (≥2 independent transitions, pre-registered pole contrasts)
  + Pf8 rare-variant boundary + official-OrthoMCL sensitivity holding D1 concordance.
- →WEAK: official OGs dissolve branch concordance (Fisher n.s.) or cons-matched controls erase
  D3 separation, or AP2 enrichment collapses under expression-matched background.

## Outputs (data/derived/WP4R/M4R-D/)
M4RD_D1_branch_substitutions.tsv (24,520 events); M4RD_D1_concordance_summary.tsv;
M4RD_D1_concordance_aux.tsv; M4RD_D2_site_architecture.tsv; M4RD_D2_partition_summary.tsv;
M4RD_D3_matched_controls.tsv; M4RD_D3_adaptive_grantham_quantiles.tsv;
M4RD_D4_zhang_constraint.tsv; M4RD_D4_essentiality_note.txt; M4RD_D5_pf_features.tsv;
M4RD_D5_regulator_enrichment.tsv; M4RD_D5_stratified_{AP2,PUF_RNA,CHROM,regulator_any}.tsv;
M4RD_D5_stratified_summary.tsv; M4RD_D6_natural_experiment.tsv (9 OGs);
claim_impact.md; README.md; params.yaml; input_manifest.tsv (+sha256).
QC: LCR masked/unmasked dual-track (convention + 2 sensitivities); within-OG permutation;
leave-one-clade-out deferred to round 2 (branch-aware replaces it as phylo control);
effect sizes + CIs throughout, never p-only.
Prohibited extrapolations respected: terminal layer kept; no mechanism language;
no repeat-first/dual-axes wording; no AlphaFold/DL/hyperparam; M4R-X01 untouched.
