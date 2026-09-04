# Gate B Review — M2 Cross-Layer Counterfactual (Codon + AA)

- gate_id: GATE_B
- review date: 2026-09-03
- decision candidates: GO / PIVOT / STOP-DOWNGRADE / HOLD
- pre-specified criterion: Gate B — stable residuals justify TRANSLATION_PRIMARY or PLASTICITY_PRIMARY as primary axis; residual robust to phylogeny, annotation, LCR and leave-one-clade-out (03 §4, 04 §9)
- supporting evidence IDs: EVID-M2-001 (M2-01 codon counterfactual, 16 spp, AA-fixed, 200 pooled perms, LCR dual-track), EVID-M2-002 (M2-02 AA counterfactual, CLR OLS genome_gc, 3647 OGs, structural/LCR dual-track)
- contradicting evidence IDs: none blocking; genome-wide adaptive codon/AA re-optimization NOT supported as compositional-maintenance DOWNGRADE on M2 alone (see §Gate)
- data quality limitations: SP009/SP010 ovale have no cds/protein (has_cds=false, frozen, n=16 for M2); codon n_perm 200 pooled-multinomial vs spec 1000 (pooled null over-dispersed vs per-gene conditional → conservative p, params.yaml N_PERM switch, seed 42, +1/BH); codon joint 6-fold null confounded by position-1 composition for L/S/R — post-Gate-B split-null audit (`M2-01_subfamily_split_check.tsv`): Leu prop_explained 0.391→0.79/0.74 within subfamilies, L downgraded to provisional; dinucleotide secondary null does not improve (median 0.397 vs 0.845, worse 16/18) — GC3-only retained; AA orthogroups fallback_python_kmerRBH_v1 (kmerRBH without diamond/MCL, USE_ALIGN=False, primary inference OG-agnostic CLR); AA site-level `PLACEHOLDER_QC_PASS` — real QC gap 1.2% ident 84% from pairwise alignments, per-site numerics empty (audit fix: synthetic placeholder values removed), full MSA/ancestor deferred pending MAFFT; AA structural stratification descriptive counts only (no cross-stratum control); N-terminal 30-codon p pooled via full; phylo GLS n=16 with distance-saturated k-mer tree (V≈diagonal ⇒ GLS≈OLS, p equal to ~6 digits — near-vacuous correction, disclosed; independence rests on LOO refits)
- sensitivity results: Codon 7/7 LOO retain gradient (median prop_explained 0.756–0.939, full 0.845; LOO reuses full-data q on subsets, significance leg near-vacuous — direction leg carries discrimination), 7 solid robust (F,I,P,A,N,D,E, ≥5/7 + direction) + L provisional (split-null audit); AA 7/7 LOO refits retain gradient (median R2 0.746–0.837, full 0.787; per-clade sig counts 14–16), 13/14 GC-coupled AAs robust (A,R,N,D,Q,G,I,K,F,P,W,Y,V + 1 sensitive S); LCR masked vs unmasked identical both layers (codon delta 0.00005, AA delta 0.00078, significance sets identical); phylo GLS 6/18 codon p<0.05 and 13/13 robust AA sign-consistent OLS/GLS (but GLS≈OLS on saturated tree — decorative, disclosed); annotation spotcheck codon 0.998 identity (150 genes SP001/005/011 NCBI vs PlasmoDB)
- independent lineages/replicates: 7 clades × dual-track; low-GC Laverania+avian and high-GC vivax+piroplasm contrasts both survive M2 controls — Gate A grade A retained
- unresolved alternatives: Codon residual GC-correlated, not expression/tRNA-correlated (M2 by design); AA residual GC-slope split (GC-rich A,R,G,P,W,V positive vs AT-rich N,D,I,K,F,Y negative) — demand axis requires M3 D_{c,s,t} + tRNA/ribo/mRNA second layer; locus-specific adaptation not excluded but not detected as genome-wide program; plasticity (indel/CNV/repeat) not tested in M2, remains fallback
- minimum next information gain task (only for HOLD): n/a — GO
- cost of next task: M3 processed-first (05/06): D006–D013 matrices/metadata small-to-moderate, no raw MS/SRA bulk by default; B-level download now justified by Gate B GO
- expert input: deep-dive second species priority (Pv P01 vs Pk H vs Pb ANKA) still open for M3 interpretation
- final decision: **GO — narrowed TRANSLATION_PRIMARY, conditional on M3**
- active axis after decision: **TRANSLATION_PRIMARY (narrowed: 7 solid codon families F,I,P,A,N,D,E + L provisional [6-fold split-null audit] + 13 GC-coupled AA families A,R,N,D,Q,G,I,K,F,P,W,Y,V as tiered candidates, not genome-wide claim)**
- claims upgraded/downgraded: CLM02 HYPOTHESIS → SUPPORTED_CANDIDATE (narrowed tiered, composition explains 84.5% median chi²-deviation at codon layer, 15% small fragmented residual, 7 solid + L provisional robust families, genome-wide adaptive codon claim DOWNGRADED without M3 second layer per 04 §7); CLM03 HYPOTHESIS → SUPPORTED_CANDIDATE (narrowed tiered, composition explains 78.7% median AA CLR R2, 13 GC-coupled robust + 6 uncoupled + 1 sensitive + 2 borderline, genome-wide adaptive AA claim DOWNGRADED without M3); CLM04 remains HYPOTHESIS (translation stage-program gated on M3 D_{c,s,t} + ≥2 orthogonal layers); CLM07 remains FALLBACK (activated only if M3 translation fails)
- project_state updated: yes — M2_COMPLETE, Gate B GO 2026-09-03, WP2/M2 closed

## Evidence Summary (WP2)

### M2-01 Codon Counterfactual (16 spp CDS, AA-fixed synonymous null)

- Background: per-gene GC3 expectation (pooled within AA family), 200 pooled-multinomial perms (multinomial N, pooled p=E/N, seed 42, +1/BH FDR per track), LCR masked/unmasked dual-track, N-term 30-codon stratum, 7 LOO, phylo GLS (V Brownian from M1-02 k21 tree). Metric `1 - chi2_GC3/chi2_uniform` per family × species.
- Result: 288 tests (16×18 families), median prop_explained **84.5% masked** (0.8451; unmasked 0.8452), 284/288 (98.6%) q<0.05 masked (power from large N, CramerV small) — 15% unexplained is fragmented. **7 solid robust** F,I,P,A,N,D,E (≥5/7 LOO + 75% direction + sig) + **L provisional** (post-Gate-B split-null audit: joint 6-fold null confounded by position-1 composition, CTN share vs genome GC r=0.990; Leu prop_explained 0.391→0.79/0.74 within subfamilies, `M2-01_subfamily_split_check.tsv`); 10 fail direction/LOO. Phylo GLS 6/18 p<0.05 GC-correlated residual (GLS≈OLS, saturated tree, disclosed). LCR delta <0.001. Spotcheck 0.998 identity (150 genes). Dinucleotide secondary null worse (median 0.397) — GC3-only retained. Tables: `M2-01_background_explained.tsv`, `M2-01_robustness_by_aa.tsv`, `M2-01_phylogeny_gls.tsv`, `M2-01_residual_table.tsv`, `M2-01_qc.tsv`, `M2-01_subfamily_split_check.tsv`, figures in `figures/`.
- Verdict: Compositional null confirmed and quantified at codon layer; genome-wide adaptive codon re-optimization not supported without expression/tRNA/mRNA second layer (04 §9 steps 1–2 pass for 8 families, steps 3+ deferred to M3). See `data/derived/WP2/M2-01_codon/claim_impact.md`, `params.yaml` spec_vs_actual_note.

### M2-02 AA Counterfactual (16 spp protein, CLR compositional)

- Background: CLR `log(p/sqrt(gmean))` + OLS `CLR ~ genome_gc` per AA (n=16, BH FDR per track), LCR 0.184% residues dual-track, structural stratification descriptive-only (KD window19>1.6, signal N30 window8>1.8, length quartiles — no cross-stratum control), orthogroups fallback_python_kmerRBH_v1 (3647 OGs ≥4 spp, 181 single-copy core, SP001-anchored RBH, kmer cosine top5, USE_ALIGN=False, REUSE_ORTHO), site layer centre-star (100 OGs sampled, QC gap 1.2% ident 84% from real pairwise alignments; per-site rows PLACEHOLDER_QC_PASS with empty numerics, parsimony ancestor DEFERRED), 7 LOO refits, phylo GLS (≈OLS, disclosed).
- Result: 20 AAs, median R2 **78.7% masked** (0.7874; mean 59.5%; unmasked 0.7866; cds_gc median 0.843 consistent), 14/20 q<0.05 masked (= unmasked) — significance of the GC-coupling slope itself (background), not residual-outlier tests. **13/20 GC-coupled robust** A,R,N,D,Q,G,I,K,F,P,W,Y,V (93% of coupled); 1 sensitive S (4/7 LOO); 6 uncoupled C,E,H,L,M,T (+ M/T borderline p≈0.055). All 20 direction-consistent masked/unmasked. LOO median R2 0.746–0.837 (coccidian-removed most sensitive -0.041). Per-species residual outliers (Cohen's d) computed but unused until M3 demand data. Structural site/OG-agnostic primary inference preserved. Tables: `M2-02_background_explained.tsv`, `M2-02_robustness_by_aa.tsv`, `M2-02_phylogeny_gls.tsv`, `M2-02_structural_stratification.tsv`, `M2-02_site_level_residual.tsv`, orthogroups cache.
- Verdict: Compositional maintenance quantified at AA layer; GC-rich vs AT-rich slope split remains compositional axis, not yet demand axis; genome-wide adaptive AA claim DOWNGRADED without M3 D_{c,s,t}+tRNA/proteome/population second layer. See `data/derived/WP2/M2-02_aa/claim_impact.md`.

### Gate B Verdict

**GO with narrowed TRANSLATION_PRIMARY candidate set, conditional on M3.** Active axis is tiered candidates (8 codon + 13 AA families), not a genome-wide adaptation claim. Composition explains ~4/5 of variation at both layers — a positive confirmation of Gate A C1 as quantitative C2. LCR is coupled but not causal. Independence (grade A) survives. Plasticity remains fallback (CLM07) if M3 translation evidence fails. M3 is required and justified (processed matrices/metadata first, ≥2 orthogonal functional layers per 02 §4 C3).

This gates B-level downloads D006–D025 per 05/06.

## What Was NOT Done (per 05/06 gating, not a deviation)

- No D006–D025 fetched at Gate B (all B/C-level: GEO/PRIDE/MCA/ENA/MalariaGEN, processed_first). D001–D005 already DOWNLOADED_VERIFIED satisfies M0–M2 A-level. Pausing B/C until Gate B GO is compliant — M3 will fetch processed matrices/metadata first, raw only if needed.
- No AlphaFold/structure ML, full-genome translation modeling, or host-omics integration (M3/M4 scope).
- No ML ranking or target nomination (mechanism precedes target per 04).

## Next

WP3/M3 supply/demand: D_{c,s,t} expression-weighted codon/AA demand + tRNA abundance/charging/modification + ribo/mRNA stability as second orthogonal layer; test whether 7(+L)+13 GC-coupled sets predict demand beyond genome_gc, cross-lineage stable; then Gate C axis lock or plasticity pivot.

## Pointers

- M2-01: `data/derived/WP2/M2-01_codon/` (run_codon_counterfactual.py, background/robustness/qc/phylo/spotcheck/figures, params, claim_impact, checksums, summary)
- M2-02: `data/derived/WP2/M2-02_aa/` (run_m2_02_aa.py, background/robustness/qc/phylo/structural/site/orthogroups, params, claim_impact, checksums, summary)
- Summary: `data/derived/WP2/GATE_B_EVIDENCE_SUMMARY.md` (Figure 2 candidate: codon residual log2 obs/exp + AA CLR slope)
- Frozen panel: `data/metadata/frozen_assembly_candidates.tsv` (18 ACC, 16 with CDS)
- Data registry: `docs/manifests/data_registry.tsv` (D001–D005 VERIFIED; D006–D025 NOT_STARTED until M3)
