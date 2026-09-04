# L1 claim impact (exploratory, D-018/D-019 — no preregistration)

- date: 2026-09-04
- verdict: **DOUBLE_POSITIVE**
- discovery hits (FULL 61-codon scan + focus extras; >=2 Dd2 backgrounds q<0.05 same sign): ['AAT_freq', 'AGT_freq', 'ATT_freq', 'CAC_freq', 'GAA_freq', 'GAG_freq', 'GAT_freq', 'GCC_freq', 'GCT_freq', 'GTC_freq', 'GTT_freq', 'TTA_freq']
- replicated in GSE226632 dTE (same sign, q<0.10): ['AAT_freq', 'GTT_freq', 'GTC_freq', 'GCT_freq', 'GCC_freq']
- rule: double-positive -> evidence card, claimable in main text; discovery-only -> supplementary, labeled single-dataset; double-negative -> H3 gene-level no-evidence.
- CLM04 impact: none until user review; Gate C PIVOT untouched.
- caveats: microarray log2ratio vs 3D7 pool; Cam3II DHA rep1-only by design (descriptive only); gc3 covariate included (beyond-composition signal); GENE-CLUSTER bootstrap B=1000 (all rows of a gene resampled together); exploratory label per D-019; ATT_freq discovery sign (negative under DHA) FLIPS in starvation dTE (positive) -> cross-perturbation sign-inconsistent, hence discovery-only.
