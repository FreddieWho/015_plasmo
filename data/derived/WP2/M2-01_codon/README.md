# M2-01 Codon Counterfactual (WP2)

**Question** In fixed-AA, fixed-length, fixed-position counterfactual, how much of synonymous codon usage variation is explained by local composition (per-gene GC3 ± dinucleotide), and what residual (direction, significance, phylogeny, robustness) remains?

**Design (04 §4, 02, 03 gate)** One primary inference level. Representative transcript = per gene locus_tag longest CDS (ties earliest); filters: length<150, not divisible by 3, internal stop, partial/frameshift excluded (max 1 CDS/gene). Dual-track LCR masked/unmasked (protein window 30 aa, Shannon entropy <1.5, `numba.njit` accelerated; codon masked if its AA lies in any LCR window). Gate B verdict uses **masked**; unmasked reported in parallel for sensitivity. 18 frozen accessions respected — SP009/SP010 `has_cds=false` excluded from CDS but retained in tree/robustness denominator. No D006–D025 / tRNA / proteome used.

**Inputs (frozen, sha in input_manifest.tsv 35 lines)**
- `cds_from_genomic.fna` 16 species (SP009/010 absent) + 2 unfixed genomic.fna for background only
- `M1-01_composition_table.tsv / M1-01_dinucleotide.tsv` (genome/GC/GC3/dinuc 16 freqs)
- `M1-02_species_tree.nwk` (Brownian V: shared root-to-MRCA branch length) + 7 clades [Laverania, vivax_knowlesi, malariae_ovale, rodent, avian, piroplasm, coccidian]
- `frozen_assembly_candidates.tsv`
- PlasmoDB-71 transcripts (SP001/005/011) for 50-gene NCBI vs PlasmoDB spot-check

**Null backgrounds (both explicit per params.yaml)**
- Primary: per-gene GC3 binomial at third position. For each gene `g` with GC3=`gc3_g`: `p(codon|AA_g) ∝ gc3_g` if codon ends G/C else `1-gc3_g`, normalized per AA family (6-fold Leu/Arg/Ser as single 6-codon family). `E_gc3 = Σ_g n_AA_g * p(codon|gc3_g)` summed per species. Stratified full / N-term 30 / core(31..end); locus fixed by construction. **Post-hoc audit (2026-09-03): joint 6-fold null confounds position-1 composition (CTN fraction of Leu vs genome GC r=0.990). Split-subfamily re-test (`run_subfamily_split_check.py` → `M2-01_subfamily_split_check.tsv`): TTR/CTN, AGY/TCN, AGR/CGN separately — Leu joint prop_explained 0.391 → 0.79/0.74 within subfamilies, so ~half of Leu's apparent residual was subfamily-shift artifact; small significant residuals (|GC-end resid| ~1.6–2.8 pp) remain. Leu candidacy downgraded to provisional.**
- Secondary: `GC3 × dinucleotide` — same `p_GC3` multiplied by per-species dinucleotide odds for codon (product of two dinuc freqs spanning codon positions / uniform 1/16) then renormalized per AA family. Captures neighbor bias beyond single-base GC. **Post-hoc audit: dinuc correction does NOT improve the null — median prop_explained 0.397 vs 0.845 for GC3-only, worse for 16/18 families (overshoots in AT genomes, e.g. AAA). GC3-only retained as primary; dinuc expected columns are descriptive only and were never used in any test.**
- Descriptive only: RSCU/ENC/PR2 computed per species for QC (ENC Wright 1990, PR2 AT/GC bias), never as adaptive evidence.

**Statistics**
- Per AA family (18 degenerate families `k>1`) chi-square `(obs-exp)^2/exp` vs GC3 null, `df=k-1`; CramerV `sqrt(chi2/(N*(k-1)))`; per-codon `residual=obs-exp`, `residual_z`, `log2(obs/exp)`; composition handled via multinomial/chi-square within family (compositional sum=1), not raw proportion Pearson.
- Permutation p: **200 pooled-multinomial draws** (spec 1000 → 200, see params `spec_vs_actual_note`). Pooled `p = E_GC3/N` per AA per species, `sim ~ Multinomial(N, pooled_p)` for `n_perm=200`, `p=(#null>=obs+1)/(n+1)`, seed 42. Heterogeneity across genes collapsed — conservative; `N_PERM` one-line switch to 1000. BH FDR per track per stratum, `q<0.05` sig. Effect + FDR + direction (±) reported together (04 §8).
- Proportion explained `1 - chi2_GC3 / chi2_uniform` per family, median across 16×18 tests.

**Phylogenetic & robustness**
- Phylo GLS: trait = per-species GC-ending residual proportion (`obs_GC_end_frac - exp_GC_end_frac`) vs predictor genome GC (alternative GC3); Brownian `V_ij` from M1-02 tree, `beta_GLS=(X'V^-1X)^-1 X'V^-1 y`, `n=16` species with CDS, reported per AA.
- LOO: leave-one-clade-out 7 clades (Laverania 3, vivax_knowlesi 4, malariae_ovale 1 with CDS, rodent 3, avian 1, piroplasm 3, coccidian 1). Metrics: median prop explained, `n_sig_masked_q05`, direction. Robust = significant `q<0.05` masked and same GC-end residual sign in ≥5/7 LOO.
- Cross-lineage replication and direction consistency ≥75% of independent clade representatives.

**Outputs**
- `M2-01_codon_residual_table.tsv` 944 rows (16×18×~ codons) — per codon observed/expected_GC3/dinuc, residual/z/log2, p/q (both tracks)
- `M2-01_background_explained.tsv` 288 rows (16×18) — chi2_GC3/uniform, prop_explained, CramerV, N, p/q per AA (masked+unmasked)
- `M2-01_gene_level_residual.tsv` 82376 genes — per-gene GC3, LCR frac, chi2, length
- `M2-01_qc.tsv` 16 rows — n_genes, retained/filter counts, codons, LCR codons, mean/median GC3, ENC/PR2 (both tracks), genome GC
- `M2-01_phylogeny_gls.tsv` 18 rows — beta/se/p_GLS, ols r/p, mean residual per AA
- `M2-01_robustness.tsv` 7 rows LOO + `M2-01_robustness_by_aa.tsv` 18 rows per-AA robustness
- `M2-01_ncbi_vs_plasmodb_spotcheck.tsv` 150 rows (50 genes × 3 species) — NCBI vs PlasmoDB identity
- `M2-01_summary.json` + `figures/*.png` (background_explained, phylo_gls, qc_enc, robustness)
- `M2-01_subfamily_split_check.tsv` (post-Gate-B audit: L/S/R split-null re-test) + `run_subfamily_split_check.py`
- `params.yaml`, `input_manifest.tsv` (35 sha), `checksums.sha256` (15 outputs), `run_codon_counterfactual.py` (49KB, numba)

**Key numbers (masked Gate B track)**
- Median prop explained by GC3: **0.845** (unmasked 0.845) — local GC is dominant null, a positive result not a negative. (Metric = per-family share of chi² deviation-from-uniform explained by GC3, not variance explained; distribution heavy-tailed: p25=0.33, 45/288 negative where GC3 underperforms uniform.)
- Residual tests: **284/288 = 98.6%** significant `q<0.05` (unmasked 283/288) — small but systematic residual after GC.
- Robust (≥5/7 LOO + direction): **8/18 AAs** — F,L,I,P,A,N,D,E. Others not robust (single-species/LCR-sensitive). **Post-hoc audit: L provisional** — joint 6-fold null misses position-1 composition (CTN share of Leu vs genome GC r=0.990); split-subfamily re-test (`M2-01_subfamily_split_check.tsv`) raises Leu prop_explained 0.391→0.79/0.74, so ~half of Leu's residual was subfamily-shift artifact; small significant within-subfamily residuals remain (~2 pp). Solid robust set = 7 (F,I,P,A,N,D,E); Leu pending M3 demand data. Caveat: the 5/7 LOO significance leg is near-vacuous (7/7 for all 18 families because 98.6% of tests are significant); effective discrimination comes from direction consistency, and codon LOO reuses full-data q-values on subsets (M2-02 refits per subset — asymmetry documented).
- Phylo GLS significant `p<0.05`: **6/18** AAs — residual correlates with genome GC after phylogeny correction; not genome-wide translation adaptation.
- Spot-check mean identity NCBI vs PlasmoDB: **0.998** (50 genes each SP001/005/011, CDS substring of transcript), annotation choice not driving signal.
- LOO range: median prop 0.756–0.939 depending on clade removed; gradient retained in all 7 runs.

**Reproduce**
```bash
LD_LIBRARY_PATH=/opt/anaconda3/lib:$LD_LIBRARY_PATH python -u data/derived/WP2/M2-01_codon/run_codon_counterfactual.py
# ~101s at n_perm=200; set N_PERM=1000 in script for spec-exact run
sha256sum -c data/derived/WP2/M2-01_codon/checksums.sha256
```

**Boundaries & non-drift**
- 18 accession invariant honored; SP009/010 explicitly `has_cds=false`.
- LCR coupling mechanistically expected (M1-03 r≈-0.82) — dual-track, Gate B masked.
- RSCU/ENC/PR2 already computed but NOT used as adaptive evidence.
- No claim beyond codon residual — amino-acid and supply/demand (D, MCA) are separate modules (M2-02 etc.). Candidates require second evidence layer per 04 §9.

**Env** python 3.11.5, numpy 1.26.4, scipy 1.13.1, pandas 2.3.3, biopython 1.87, matplotlib 3.8.4, numba (LCR). Run date 2026-09-03, seed 42.
