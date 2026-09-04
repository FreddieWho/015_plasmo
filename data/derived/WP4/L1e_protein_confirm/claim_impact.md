# L1e claim impact (ring 3: protein layer; exploratory)

- date: 2026-09-04
- data: Small-Saunders 2024 (Nat Microbiol) TMT MSstats source data, Dd2 WT vs R539T, 12h post DHA/DMSO pulse.
- prediction from chain (L1 mRNA up + GSE226632 TE up + Asn supply down): Asn-rich proteins shift under DHA.
- results per comparison in L1e_protein_tests.tsv; chain verdict updated in STATUS after user review.
- caveat: TMT detects ~2-3k proteins (coverage biased to abundant); Asn-rich proteins are often LOW abundance
  (regulatory) -> detection bias against the test set, noted for honest interpretation.
- claim impact: protein-layer evidence for/against the Asn-axis chain; no claim change yet.


## Verdict (2026-09-04): NOT CONFIRMED at protein layer with this dataset
- Direct test WT_DHA vs WT_DMSO @12h: Asn-top10% vs rest MW p=0.28 (n.s.); spearman rho=+0.069 (p=0.014, tiny).
- Strongest signal is a GENOTYPE contrast (WT vs R539T under DHA: Asn-rich proteins more abundant in WT, p~0)
  - steady-state genotype difference, not a stress response.
- Structural limitations: TMT covers only 1249/5274 genes (24%); Asn-rich regulators nearly invisible
  (107/528 Asn-top10 detected; AP2: 2/28; PUF: 0). Detection bias + 12h integration time vs <=6h mRNA signal.
- Consequence: chain closure NOT claimed; L-007 stays open pending starvation-resolved proteomics (Mok 2021).
- This negative is reported, not rescued (per adversarial-review rule 7).
