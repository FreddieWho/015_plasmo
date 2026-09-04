# L1: gene-level synonymous-codon x stress interaction (exploratory)

- direction: D-018 / D-019 (exploratory; no preregistration; replication gate = only path to claim)
- date: 2026-09-04 | package: PLASMODIUM-C2F, post-Gate-C (PIVOT) direction L1

## What
Tests whether a gene's synonymous codon composition predicts its stress response at the GENE level (the statistically correct resolution; M3-03's aggregate demand-share test was insensitive by construction).

- Discovery: GSE151189 (Mok 2021) — per-gene DHA response per K13 background (Dd2 WT/R539T/C580Y), early times (<=6h), model `resp ~ feature + gc3 + time`, bootstrap B=1000 over genes.
- Replication: GSE226632 (Li 2024) — per-gene TE change under amino-acid starvation (dTE = TE_AA_free - TE_CM), model `dTE ~ feature + gc3`.
- Focus features (mechanism-anchored): AAA_freq, AAG_freq, AAA_share, ATT_freq, ATCATA_freq; Lys pair from Small-Saunders 2024 s2U@U34 mechanism, Ile from M2 robust families.
- Exploration: all 61 sense-codon frequencies, labeled exploratory.
- LCR dual-track for focus features (unmasked + masked).

## Verdict rule
double-positive (discovery q<0.05 same sign in >=2 Dd2 bgs AND replication same sign q<0.10) -> claimable; discovery-only -> supplementary; double-negative -> H3 gene-level no-evidence.

## Key caveats
- microarray log2 ratios vs 3D7 pool (relative, not TPM)
- Cam3II DHA arms are rep1-only at every timepoint (experimental design) -> Dd2-only inference; Cam3II descriptive
- gc3 covariate included throughout: signal must be BEYOND composition background
- POST-hoc provenance (follows M3-03 null) -> exploratory label per D-019

## Files
- L1_response_by_gene.tsv — per-gene DHA response per bg x time
- L1_focus_slopes.tsv / L1_exploration_slopes.tsv — bootstrap slopes + BH q (per track)
- L1_GSE226632_dTE.tsv / L1_TE_replication.tsv — replication layer
- L1_K13_row.tsv — K13 (PF3D7_1343700) feature values + responses (control)
- L1_qc.json, claim_impact.md, params.yaml, input_manifest.tsv, checksums.sha256
