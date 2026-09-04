# Gate C Review — M3 Functional & Perturbation Integration

- gate_id: GATE_C
- review date: 2026-09-04
- decision candidates: GO / PIVOT / STOP-DOWNGRADE / HOLD
- pre-specified criterion: Gate C — at least one module with >=2 orthogonal support layers and core phenotype link (03 §5); D-016 rule: M3-03 early DHAxK13 hit(q<0.05, LysAAA/AAG+Ile, direction-consistent) -> GO-narrowed-Pf, else PIVOT (user pre-approved 保A with auto-fallback)
- supporting evidence IDs: EVID-M3-001 (M3-01 demand/TE: Pf-internal 271/405 q<0.05, split-subfamily gains R_AGR+0.56/L_CTN+0.20), EVID-M3-002 (M3-02: s2U/LysAAA-AAG/K13/RSA-drug panel, GSE151189 VERIFIED), EVID-M3-003 (M3-03: assay-sensitive null — K13 mRNA +0.4..+1.0 log2 under DHA proves detection power)
- contradicting evidence IDs: EVID-M3-001 primary Pf-vs-Pk 13/36 p=0.13 (cross-species program fails); EVID-M3-003 primary 0 early hits + 0/8 TE-sign consistency (DHA-perturbation contrast fails at demand layer)
- data quality limitations: GSE151189 microarray two-color log2ratio vs 3D7 pool (relative proxy 2^w, not absolute TPM); DHA n=3v3 per Dd2 stratum (large-effect power only), Cam3II 8h DHA n=1 (NA, excluded from primary); late-time stage-delay confounding (early<=8h primary pre-registered); tRNA charging file empty (regex gap, disclosed); TMT quants still in pending xlsx; HS_ tRNA human-annotated proxy
- sensitivity results: BH per bg×time; effect sizes negligible (demand-share deltas ~1e-4); best raw p Dd2_WT 3h ATT 0.0062 (q=0.13), DiD ATC/ATT 0.051/0.056 (q=0.19); K13 mRNA UP under DHA (+0.4..+1.0) vs K13 protein DOWN per M3-02 => regulation is post-transcriptional (supports translation-layer framing, not demand-shift framing)
- independent lineages/replicates: second-lineage requirement NOT met (Pk fails primary, Pb aux descriptive only, Pv population-only); strong-perturbation contrast (DHA×K13, 5 backgrounds × time × 3rep) does NOT rescue: 0 hits
- unresolved alternatives: tiny effects below detection (1e-4 shares) cannot be excluded — but biologically negligible for a "program" claim; stage-label mismatch across species remains as alternative for M3-01 primary (moot after M3-03, which is within-Pf)
- minimum next information gain task (only for HOLD): none — HOLD already consumed (M3-03 was the HOLD task per D-013)
- cost of next task: WP4 narrowed Pf mechanism focus uses in-hand data (M3-01/02/03 + D008 PDFs); no new downloads; T01 journal decision pending user
- expert input: requested for Pf-specific K13/s2U mechanism novelty + validation feasibility (Gate D input moved forward)
- final decision: **PIVOT — species-specific (Pf-centered LysAAA/s2U/K13 translation mechanism); cross-species program claim DOWNGRADED**
- active axis after decision: **TRANSLATION_PRIMARY, Pf-specific only**
- claims upgraded/downgraded: CLM04 HYPOTHESIS → SUPPORTED_CANDIDATE narrowed Pf-only (within-Pf TE/decay + tRNA-mod + drug panel + K13 mRNA = >=2 orthogonal layers + ART core phenotype; cross-species wording DOWNGRADED); CLM02/CLM03 unchanged (compositional maintenance stands, Figure 2 as calibration per D-015); CLM07 stays FALLBACK (not activated — translation axis survives in Pf-specific form); CLM05/CLM06 remain HYPOTHESIS (gated on M4/Gate D)
- project_state updated: yes — M3_COMPLETE_GATE_C_PIVOT, M4/WP4/GATE_D next, T01 journal re-assessment pending user approval

## Evidence Summary (WP3)

### M3-01 Demand (EVID-M3-001)
Pf-internal TE/decay strongly codon/AA-associated (271/405 q<0.05, rep corr 0.92-0.98); split-subfamily null gains; PRIMARY Pf-vs-Pk 13/36 p=0.13 NEGATIVE. Post-hoc residual-vector correlation codon 0.88/0.43 (large-point driven), AA -0.59/-0.52 (mirror = composition axis). Tables + 04-s10 deliverables in `data/derived/WP3/M3-01_demand/`.

### M3-02 Functional extraction (EVID-M3-002)
mcm5s2U+mCm down in ART-R+DHA; s2U@U34 Lys/Glu/Gln; PfMnmA cKD 4%->11% survival; LysAAA/AAG top driver, K13 52/57 AAA; drug panel AZT/FSM 2x sens, DSM265 unchg, ATQ plateau, 42C unchg. GSE151189 VERIFIED, PXD043747 all-RAW no-download. In `data/derived/WP3/M3-02_functional/`.

### M3-03 DHA x K13 (EVID-M3-003, the HOLD task)
156 samples, 5 backgrounds, probe->PF3D7 direct (10,214 probes, 5,107 genes). Primary: 0 early hits in LysAAA/AAG+Ile set; DiD all q n.s.; demand deltas ~1e-4 (negligible); TE-sign consistency 0/8. K13 mRNA UP +0.4..+1.0 (assay-sensitive null: detection works, biology doesn't move at demand layer). In `data/derived/WP3/M3-03_dha/`.

## What WAS Done (gating compliance)
- No raw SRA/MS bulk at any point (largest fetch: MCA 21MB + GEO 9.7MB, all processed)
- No ML ranking, no target list, no host omics, no wet-lab proposal
- LCR dual-track throughout; phylogeny-respecting units; negatives recorded, no subgroup dredging (Cam3II n=1 strata excluded, not imputed)

## Next
WP4/M4 narrowed Pf mechanism: 1 main mechanism (AAA/s2U/K13 translational persistence) + 1-3 nodes; two competing models (transcriptional vs translational control of K13 output); essentiality/selectivity/escape bounds from public data; T01 journal decision (B-track: Nat Commun / PLoS Pathog / iScience) awaits user.

## Pointers
- M3-01: `data/derived/WP3/M3-01_demand/`; M3-02: `data/derived/WP3/M3-02_functional/`; M3-03: `data/derived/WP3/M3-03_dha/`
- Data registry: D006/D010/D011/D012/D013 VERIFIED, D008 PARTIAL, D009 meta-only
- Gate C criterion source: docs/03_ROADMAP_GATES_AND_DELIVERABLES.md §5; rule D-016
