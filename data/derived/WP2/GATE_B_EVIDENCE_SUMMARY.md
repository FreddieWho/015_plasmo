# Gate B Evidence Summary — WP2 Cross-Layer Counterfactual (M2-01 + M2-02)

- date: 2026-09-03
- gate: GATE_B (M2 → axis selection)
- inputs: M2-01 codon (16 spp CDS, per-gene GC3/dinucleotide null, 200 pooled-multinomial perms, LCR dual-track, 7 LOO, NCBI vs PlasmoDB spotcheck) + M2-02 AA (CLR, OLS genome_gc, 3647 OGs 181 core, 7 LOO, phylo GLS, 0.18% LCR)
- tree: M1-02 k21 crc32 stride20 NJ rooted Tg (frozen invariant)

## 1. Headline numbers

| Layer | Background explained (median, masked) | GC-coupling significance (q<0.05 masked) | Robust (≥5/7 LOO + direction-consistent) | LCR masked vs unmasked delta | Phylo-corrected signal |
|---|---|---|---|---|---|
| **Codon (synonymous, AA-fixed)** | 84.5% per family GC3 (`1 - chi2_GC3/chi2_uniform`, 288 tests; share-of-chi² metric, not variance explained) | 284/288 (98.6%) — fragmented, small CramerV | **7 solid** F,I,P,A,N,D,E + **L provisional** (post-Gate-B split-null audit: joint 6-fold null confounded by position-1 composition, CTN share vs GC r=0.990; Leu prop_explained 0.391→0.79/0.74 within subfamilies) | 0.8451 vs 0.8452 (identical) | 6/18 phylo GLS p<0.05, GC-correlated residual; GLS≈OLS (saturated k-mer tree, V≈diagonal — disclosed) |
| **AA (CLR, compositional)** | 78.7% `R2 CLR~genome_gc` (20 AAs, mean 59.5%) | 14/20 (70%) — this is significance of the GC coupling (slope), i.e. the background itself | **13 GC-coupled robust** A,R,N,D,Q,G,I,K,F,P,W,Y,V + 1 sensitive S (per-species residual outliers computed but unused until M3) | 0.7874 vs 0.7866 (delta 0.0008) | Sign-consistent OLS/GLS for 13 robust; p_gls == OLS p to ~6 digits (near-vacuous correction, disclosed) |

Both layers: gradient retained in 7/7 LOO (codon median 0.756–0.939; AA median 0.746–0.837). N50/Tg artefact not driving. Annotation spotcheck codon 0.998 identity (150 genes SP001/005/011).

## 2. What is confirmed

- **Compositional null is the dominant process.** Single predictor (GC3 or genome GC) explains ~4/5 of variation at both layers, survives AA-fixed / CLR closure-safe / LCR masking / structural stratification / phylogeny GLS / 7 LOO. This is a positive confirmation, not a null-failure. Extends Gate A C1 to quantitative C2.
- **LCR is coupled but not causal.** M1-03 r≈-0.82 already predicted; dual-track delta <0.001 and significance sets identical. Gate B uses masked per contract, but headline unchanged.
- **Independence holds.** Low-GC (Laverania + avian) and high-GC (vivax clade + piroplasm) remain independent after M2 controls.

## 3. What remains as candidate (fragmented)

- Codon: 15% unexplained chi-square fragmented across families; 7/18 solidly robust (F,I,P,A,N,D,E) + L provisional after split-null audit, 10 fail direction/LOO. Not a genome-wide adaptive codon program. Codon LOO caveat: significance leg near-vacuous (all families 7/7), direction leg carries discrimination; LOO reuses full-data q on subsets (no per-subset refit).
- AA: 6/20 GC-indifferent (C/E/H/L/M/T), 1 sensitive (S), 2 borderline (M/T). 13 robust AAs are the **GC-coupled set** — GC-rich (A,R,G,P,W,V) positive slope vs AT-rich (N,D,I,K,F,Y) negative — a compositional axis, not yet a demand axis without `D_{c,s,t}`. Structural stratification computed as descriptive counts only (no cross-stratum control; earlier "survives stratification" wording corrected).
- Both signals are **GC-correlated**, not expression/tRNA-correlated (not tested in M2 by design). Cannot license C3 translation-program claim.

## 4. Translation vs plasticity at Gate B

| Axis | M2 support | M2 against | Gated data still needed |
|---|---|---|---|
| **TRANSLATION_PRIMARY** | 7 solid codon families (+L provisional) + 13 GC-coupled AA families pass 04 §9 steps 1–2 (stable, not single-clade/LCR) — valid tiered candidates for M3 `D_{c,s,t}` + tRNA/ribo/mRNA test | Genome-wide adaptive codon/AA re-optimization NOT supported on M2 alone; effect sizes small; phylo signal is GC, not demand | M3 D006–D013 supply/demand, ≥2 orthogonal functional layers (02 §4 C3) |
| **PLASTICITY_PRIMARY** | No M2 test — indel/CNV/repeat not analysed by design | No evidence to select it now | M3 NAR/ENA + MalariaGEN/EVOL indel/CNV (D015, D020–21) |
| **DOWNGRADE/STOP** | Would discard stable 7(+L)+13 candidates and stop before testing their predicted function | Not warranted — stable residuals meet GO-candidate threshold for follow-up | — |

## 5. Recommendation for Gate B (single axis per 03 §4 hard rule)

**GO with narrowed TRANSLATION_PRIMARY candidate set, conditional on M3.**

- `active_axis = TRANSLATION_PRIMARY` (narrowed: 7 solid codon families F,I,P,A,N,D,E + L provisional + 13 GC-coupled AA families A,R,N,D,Q,G,I,K,F,P,W,Y,V as tiered candidates, not genome-wide claim).
- Genome-wide adaptive codon/AA program claim is **DOWNGRADED** to compositional maintenance at C2 without M3 second layer (04 §7 prohibits L2/L3 upgrade to find signal).
- Plasticity remains fallback (CLM07) if M3 translation evidence fails.
- M3 is required and justified (cost: processed matrices/metadata first, no raw MS/SRA bulk by default).

## 6. Risks & mitigations carried to M3

- Codon `n_perm 200 pooled` (spec 1000) — pooled null is over-dispersed vs per-gene conditional null → conservative p; `N_PERM` one-line switch to 1000 (~500s), seed 42, +1/BH documented in params.yaml.
- Codon 6-fold subfamily confound (L/S/R position-1) — quantified by post-hoc split-null check (`M2-01_subfamily_split_check.tsv`); L provisional; future reruns should default to subfamily-aware nulls. Dinucleotide secondary null does not improve (median 0.397 vs 0.845) — GC3-only retained.
- AA `kmerRBH_v1` orthogroups (no diamond/MCL/MAFFT) — primary inference OG-agnostic CLR; site layer `PLACEHOLDER_QC_PASS` (real QC gap 1.2% ident 84%; per-site numerics empty, not mechanism-licensed).
- Structural stratification is descriptive only (no cross-stratum statistical control performed).
- PGLS near-vacuous at this divergence (saturated k-mer tree, V≈diagonal, GLS≈OLS); independence claims rest on LOO refits, not GLS.
- N-terminal 30-codon stratum p pooled via full — M3 can re-evaluate stratified p if codon demand localises N-terminally.
- 16 CDS species limit PGLS power — reported alongside OLS.

## 7. Outputs

- M2-01: data/derived/WP2/M2-01_codon/ (residual/background/robustness/qc/phylo_gls/spotcheck/figures + params/README/claim_impact/checksums)
- M2-02: data/derived/WP2/M2-02_aa/ (residual/background/robustness/qc/phylo_gls/structural/site/orthogroups/figures + params/README/claim_impact/checksums)
- This summary → docs/GATE_B_REVIEW.md → manifests

## 8. Figures (Figure 2 candidate)

- M2-01: residual log2(obs/exp) vs GC3 expectation + robustness by AA (7 solid + L provisional highlighted)
- M2-02: CLR R2 distribution + per-AA slope vs GC (13 robust)
