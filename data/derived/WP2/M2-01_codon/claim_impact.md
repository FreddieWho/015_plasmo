# M2-01 Claim Impact — Codon Counterfactual (Gate B input)

## Summary verdict
**Local composition explains ~84.5% (median) of synonymous codon variation.** This is a **positive compositional null confirmation**, not a failure to find signal. A systematic residual remains (98.6% of family-tests significant after GC3, `q<0.05` masked), but it is **small in magnitude, fragmented across families, and only partially robust** — does not license a genome-wide adaptive codon re-optimization claim at Gate B on codon data alone (04 §9 requires independent second layer + cross-lineage stability).

## Evidence

- **Background** `M2-01_background_explained.tsv` 288 tests (16 spp × 18 AA families):
  median `prop_explained_GC3_masked = 0.8451` (unmasked 0.8452), CramerV small. Metric = per-family share of chi² deviation-from-uniform explained by GC3 (`1 - chi2_GC3/chi2_uniform`), **not** variance explained; distribution heavy-tailed (p25=0.33, 45/288 negative where GC3 underperforms uniform). Removing any clade: LOO median 0.756–0.939 (7 runs), gradient retained in all.
- **Residual significance** (200 pooled-multinomial perms, BH FDR per track):
  `284/288 (98.6%) masked q<0.05`, `283/288 unmasked`. Significance reflects large N per family (power), not large effect — CramerV and `log2(obs/exp)` modest. 15% unexplained chi-square is fragmented.
- **Robustness** `M2-01_robustness_by_aa.tsv`:
  **ROBUST (≥5/7 LOO + direction-consistent ≥75%) = 8/18**: F, L, I, P, A, N, D, E.
  10/18 fail (direction flips or drops below significance when a clade removed). Hence per-family follow-ups possible, but no universal codon adaptation pattern.
  **Post-Gate-B audit correction (2026-09-03):** (i) the 5/7 LOO-significance leg is near-vacuous — all 18 families score 7/7 because 98.6% of tests are significant; effective discrimination comes from the direction-consistency leg alone, and codon LOO reuses full-data q-values on species subsets (no per-subset permutation refit; M2-02 does refit — asymmetry disclosed). (ii) **L downgraded to provisional**: joint 6-fold null misses position-1 composition (CTN share of Leu vs genome GC r=0.990); split-subfamily re-test (`M2-01_subfamily_split_check.tsv`, TTR/CTN separate) raises Leu prop_explained 0.391→0.79/0.74 — ~half of Leu's residual was subfamily-shift artifact. Small significant within-subfamily residuals persist (median |GC-end resid| 1.6–2.8 pp; TTR direction +15/-1), so Leu is not excluded, but its robust-set membership no longer stands on the joint-null evidence. Solid robust set = **7** (F,I,P,A,N,D,E) + L provisional.
- **Phylogeny-controlled residual** `M2-01_phylogeny_gls.tsv` (Brownian V from M1-02, n=16):
  `phylo_GLS p<0.05 in 6/18 AAs`. Residual GC-ending bias covaries with genome GC after phylogenetic correction — points to shared compositional process, not translation-driven supply/demand (which would predict expression/tRNA-correlated residual, not GC-correlated). **Disclosure:** the M1-02 k-mer tree is distance-saturated, so V ≈ diagonal and GLS ≈ OLS (p values agree to ~6 digits); "survives phylogenetic correction" is true but the correction is near-vacuous at this divergence scale.
- **LCR stratification**: masked vs unmasked nearly identical (median 0.8451 vs 0.8452; sig 284 vs 283). LCR is mechanistically GC-coupled (M1-03 r≈-0.82) but does not drive the headline number; dual-track retained and **Gate B uses masked** per contract.
- **Annotation robustness** `M2-01_ncbi_vs_plasmodb_spotcheck.tsv` (150 genes):
  mean identity **0.998** NCBI CDS vs PlasmoDB-71 transcripts for SP001/005/011 (50 genes each, CDS substring of transcript). Representative-transcript rule (longest CDS) not inflating residual.

## What is supported / weakened

| Claim thread | Impact |
|---|---|
| **CLM: genome-wide composition maintenance (Gate A)** | **Supported & quantified** — 84.5% median codon chi²-deviation-from-uniform explained by per-gene GC3 alone. **Post-hoc correction:** the dinucleotide secondary null does NOT improve on GC3 (median 0.397, worse for 16/18 families — overshoots in AT genomes); GC3-only retained. Composition gradient survives stringent within-AA, within-gene counterfactual. |
| **CLM: LCR artifact drives codon pattern** | **Weakened** — LCR masked vs unmasked identical; not an artifact driver. |
| **CLM: genome-wide adaptive codon re-optimization (e.g., tRNA matching)** | **Not supported as genome-wide inference from codon data alone.** Robust residual limited to 8/18 families and 6/18 phylo hits; effect sizes small; no expression/tRNA layer in M2-01 (by design). Upgrading to L2/L3 models to "find signal" prohibited (04 §7). |
| **CLM: specific AA families as candidates** | **Narrowed** — F,I,P,A,N,D,E pass 04 §9 steps 1–2 (residual direction stable, not single-clade/LCR driven); **L provisional** (split-null audit: ~half its residual was 6-fold subfamily-shift artifact, small within-subfamily residual remains). Require step 3+ (second evidence: expression-weighted demand, tRNA, ribo/mRNA stability) in downstream modules before mechanism entry. |

## Gate
- **Gate B on codon counterfactual alone: GO with compositional null, NO-GO for genome-wide adaptive codon claim.** The 8 robust AAs may enter tiered candidate tracking (04 §9) but must not be presented as C3-level translation adaptation without independent supply/demand evidence (M2-02, M3).
- Negative-result boundary stated: 15% unexplained is real but small and fragmented; absence of genome-wide adaptation signal here does not exclude locus-specific effects detectable with expression/tRNA integration.

## Limitations & review risk
- `n_perm 200 pooled-multinomial` approximates within-gene heterogeneity (pooled `p = E/N`, single multinomial per AA per species per perm) vs spec 1000 per-occurrence draws. Direction of bias: the pooled null is over-dispersed relative to the per-gene conditional null (Jensen: Σ n_g p_g(1-p_g) ≤ N p̄(1-p̄)), so p-values are conservative (too large). Reviewer challenge anticipated — mitigated by explicit `spec_vs_actual_note` in `params.yaml`, seed 42, `+1` smoothing, BH FDR across 288 tests, and min p = 1/201 = 0.005 (did not bind: 284/288 significant anyway).
- 6-fold families (L/S/R) share position-1 differences across subfamilies; joint null confounds this. Post-hoc split-null re-test quantified the effect (L 0.391→0.79/0.74, R 0.201→0.95/0.73, S 0.78→0.94/0.86); L provisional pending M3. Future reruns should use subfamily-aware nulls by default.
- 16 species with CDS only; SP009/010 excluded by `has_cds=false` (frozen invariant).
- N-terminal 30-codon stratification computed for background table; permutation p currently pooled via full stratum (core vs N-term signal conflated — acceptable for Gate B, noted in params).
- LOO robustness reuses full-data permutation q-values on species subsets (no per-subset refit), and its significance leg is near-vacuous (all families 7/7); direction-consistency leg carries the discrimination. M2-02 refits per subset — asymmetry disclosed.
- Phylo GLS uses distance-saturated k-mer tree (V ≈ diagonal), so GLS ≈ OLS; "phylogeny-robust" should be read as OLS-robust at this divergence scale.
- gc3 per gene is computed including the tested codons themselves (self-conditioning) and, for the masked track, including masked residues — standard approximations, negligible at these N.
- RSCU/ENC/PR2 computed but unused for inference per 04 §4.

## Downstream
M2-02 amino-acid counterfactual, then M3 supply/demand (`D_{c,s,t}`) must test whether the 8 robust families predict expression-correlated codon demand beyond GC, with tRNA/ribosome/mRNA stability as second layer.
