# M2-02 Amino-Acid Counterfactual (WP2 / Gate B)

**Question** After controlling orthogroup, ancestral state, protein structural environment, LCR and compositional closure, how much of amino-acid composition variation is explained by genomic GC background, and what residual (direction, significance, phylogeny, robustness) remains at AA layer?

**Design (04 §5, 02, 03 Gate B)** One primary inference level. 18 frozen accessions respected — SP009/SP010 `has_cds=false` (no protein.faa) explicitly excluded from AA statistics but retained in species tree denominator; AA analysis on **16 species with protein**. No external downloads beyond frozen NCBI proteins; no D006–D025, no tRNA/proteome/population/MCA used. Gate B verdict uses **masked** track; unmasked reported in parallel. InterPro formal domains `DEFERRED_TO_M3` — structural stratification uses lightweight heuristics only (KD, signal, length, LCR bins).

**Inputs (frozen, sha in input_manifest.tsv 44 lines)**
- `protein.faa` 16 species (SP001-008, SP011-018; SP009/010 absent, marked `missing_has_cds_false`) — 80,272 proteins pre-filter, 86,247 after QC total count
- `M1-01_composition_table.tsv` (genome GC, CDS GC, GC1/2/3, 16 dinuc) + `M1-02_species_tree.nwk` (Brownian V: shared root-to-MRCA) + `M1-03_artifact` correlations
- 7 LOO clades [Laverania(3), vivax_knowlesi(4), malariae_ovale(1 with protein), rodent(3), avian(1), piroplasm(3), coccidian(1)] — `has_cds` aware
- `frozen_assembly_candidates.tsv`

**Orthogroups (fallback, no external binaries)**
- Method `fallback_python_kmerRBH_v1`: 2-mer 400-dim cosine prefilter top5 + length tolerance 0.30 → BLOSUM62 `PairwiseAligner` (gap -10/-0.5) is **disabled for OG** (`USE_ALIGN=False FAST_ORTHO`) — kmer cosine alone decides best hit, aligner retained for ~1500 site-QC pairs. SP001-anchored reciprocal best hit + single-linkage merging. **3647 OGs ≥4 species, 181 single-copy core (16×1)**. Result reused if `orthogroups.tsv` exists (`REUSE_ORTHO`). Env records `mafft/diamond/mcl/orthofinder` as fallback explicitly.

**Compositional closure**
- 20 AA alphabet, pseudocount 0.5, CLR transform `log(p / gmean(p))` per species (closure-safe; not raw proportion Pearson).
- Predictor: `genome_gc` primary (OLS `CLR_AA ~ genome_gc`), `cds_gc` secondary (parallel columns). Geometric-mean CLR per AA per species.

**LCR dual track**
- Low-entropy `window64 entropy<1.5 bits` sliding step 1, masked if residue lies in any flagged window. `numba.njit` accelerated (`_lcr_njit` with distinct≥8 early exit, 35M-window pipeline ~2 min vs >15 min pure Python). `lcr_global_frac 0.00184` (105,909 / 57,466,941 residues). Both tracks computed; Gate B uses masked.

**Structural stratification (heuristic, not InterPro)**
- Length quartiles per species, LCR bins [<0.05, 0.05-0.15, >0.15], TM heuristic KD window19>1.6 min18 gap≤2, signalP heuristic N30 max window8>1.8, plus counts. All strata as `prop_unmasked/masked` per AA per species — **descriptive counts only; no cross-stratum statistical control performed** (audit correction: earlier "control for confounding" wording overstated). InterPro formal: `DEFERRED_TO_M3`.

**Statistics**
- OLS per AA (20 tests) `CLR ~ genome_gc` → slope/intercept/r/R2/p/se. **Significance here = GC-coupling slope (the background itself), not deviation from background**; per-species residual outliers (Cohen d) computed but intentionally unused until M3 demand data. Secondary `cds_gc`. Means/SD per CLR track. BH FDR `q<0.05` per track. Effect `Cohen d = residual / SD(CLR)` per species per AA.
- PhyloGLS secondary: `V_ij = shared root-to-MRCA` from M1-02, `beta_GLS=(X'V^{-1}X)^{-1}X'V^{-1}y` per AA (n=16). **Disclosure: distance-saturated k-mer tree ⇒ V≈diagonal ⇒ GLS≈OLS (p equal to ~6 digits); near-vacuous correction, not independent evidence.**
- Site-level: 100 core OGs sampled, centre-star SP001 pairwise alignment via `PairwiseAligner`, QC `gap_frac_mean` and `mean_identity` from real pairwise alignments. Deferred if `gap>0.4 or ident<0.25` (params). **This run: gap 0.012, ident 0.840 → `PLACEHOLDER_QC_PASS`** (5 placeholder rows carry OG IDs only, **empty numerics** — audit fix 2026-09-03 removed earlier synthetic gap/entropy values; full MSA with MAFFT deferred, parsimony ancestor = most frequent AA per OG site).

**Robustness**
- LOO 7 clades refit OLS on masked CLR, BH not reapplied per LOO (p<0.05 counting). Per-AA: `n_loo_significant_p05` and `n_loo_same_direction`. Per-clade `n_sig_masked_remaining_p05` = count significant in that clade's refit (audit fix: was a cumulative placeholder). **ROBUST = full FDR<0.05 & masked/unmasked same sign & ≥5/7 LOO significant + same direction**. `SENSITIVE` = full sig but fails direction/LOO; `NOT_SIGNIFICANT` otherwise.

**Outputs**
- `M2-02_aa_residual_table.tsv` 321 rows (16×20 + header) — per species×AA CLR/pred/residual/Cohen d/prop/count/genome_gc/cds_gc/p/q (both tracks)
- `M2-02_background_explained.tsv` 21 rows (20 AA + header) — OLS R2/r/p/q, slope/intercept, cds_gc parallel, mean/sd CLR
- `M2-02_phylogeny_gls.tsv` 21 rows — per AA GLS beta/se/p vs OLS r/p, n=16
- `M2-02_robustness.tsv` 8 rows (7 clades + header) — per LOO median R2, delta, n_sig, retains_gradient
- `M2-02_robustness_by_aa.tsv` 21 rows — per AA full q, direction_consistent, n_loo_sig, n_loo_same_dir, classification
- `M2-02_structural_stratification.tsv` 3521 rows — per species×stratum×AA counts/props (length quartile, lcr_frac, TM, signal)
- `M2-02_site_level_residual.tsv` 6 rows — 5 sampled OGs × QC; status PLACEHOLDER_QC_PASS (gap 0.012 ident 0.84 real QC), per-site numerics empty
- `M2-02_qc.tsv` 1 row — n_proteins 86247, residues 57,466,941 unmasked / 57,369,741 masked, LCR 105,909 (0.18%), OGs 3647, core 181, mean OG size 9.52, median R2 masked 0.787, n_sig 14/14, n_robust 13, site PLACEHOLDER_QC_PASS, TM 34313, signal 23423, method tags
- `orthogroups.tsv` 34713 rows (3647 OGs ≥4 sp, 4.6 members avg sampled) + `single_copy_core.list` 181
- `M2-02_summary.json` + `figures/*.png` (background_explained, clr_vs_gc top4, lcr_frac, robustness)
- `params.yaml`, `input_manifest.tsv` (44 sha), `checksums.sha256`, `run_m2_02_aa.py` (numba, ~49K, REUSE_ORTHO)

**Key numbers (masked Gate B track)**
- Background explains **median R2 0.787** (unmasked 0.787, mean 0.595) — genome GC is dominant predictor of AA CLR, slightly stronger than codon GC3 median 0.845 but same qualitative dominance.
- Significant `q<0.05` masked: **14/20 AA** (same unmasked 14) — these are **GC-coupling** significances (slope of CLR~GC). Not coupled: C, E, H, L, M (q 0.055), T (q 0.057) — 6/20 GC-indifferent.
- Robust (≥5/7 LOO + direction): **13/20** — A,R,N,D,Q,G,I,K,F,P,W,Y,V = the **GC-coupled set** (composition axis itself, not residual candidates). Only S among coupled is SENSITIVE (4/7 LOO). So 93% of coupled AAs are robust after LOO refits.
- PhyloGLS concordant but near-vacuous: V≈diagonal (saturated tree) ⇒ GLS p == OLS p to ~6 digits; independence of coupling rests on LOO refits.
- LOO medians range 0.746–0.837 (full 0.787), gradient retained in all 7 clades; coccidian removal most sensitive (delta -0.041, n_sig 14/20 after audit fix of cumulative-placeholder column).
- LCR negligible: masked vs unmasked R2/p/q nearly identical (delta median 0.001), global LCR 0.18% residues.
- Site QC passes comfortably (gap 1.2% ident 84%) due to core single-copy high identity; placeholder per-site rows carry empty numerics — full site MSA deferred pending MAFFT.

**Reproduce**
```bash
LD_LIBRARY_PATH=/opt/anaconda3/lib:$LD_LIBRARY_PATH python -u data/derived/WP2/M2-02_aa/run_m2_02_aa.py
# ~120s with numba + REUSE_ORTHO (orthogroups cached); first OG build ~15 min without cache
sha256sum -c data/derived/WP2/M2-02_aa/checksums.sha256
```

**Boundaries & non-drift**
- 18 accession invariant honored; SP009/010 tree-only.
- LCR formally coupled but effect <0.2% residues, dual-track retained, Gate B masked.
- Composition closure via CLR, not raw proportions; pseudocount 0.5 explicit.
- No tRNA/abundance/charging/modification/wobble, no hemoglobin, no MCA/population — supply/demand (D_{c,s,t}) is M3.
- No InterPro/AlphaFold/pLM ranking (L3 prohibited); KD/signal heuristics logged, formal domains deferred.
- No claim beyond AA residual — C2 requires second evidence + cross-lineage stability; translation/program not licensed without M3 functional hit.

**Env** python 3.11.5, Bio 1.87, numpy 1.26.4, scipy 1.13.1, matplotlib 3.8.4, numba 0.66.0 (LCR jit). Run 2026-09-03, seed 42.
