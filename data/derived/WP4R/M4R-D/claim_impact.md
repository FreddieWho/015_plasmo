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
