# L1g claim impact (ring 3b: starvation proteomics; exploratory)

- date: 2026-09-04
- data: Li 2024 iScience (PMC11544085) Table S2, NF54 6h AA starvation vs CM, 1148 proteins DE table (limma).
  SAME study as our TE replication dataset GSE226632 -> mRNA/TE/protein cross-layer test within one experiment.
- tests: Asn-top10% (genome-wide threshold) vs rest on protein logFC (MW, two-sided); spearman(Asn, logFC);
  sig-set direction fractions; stage-confound control = exclude schizont-peak genes (L1f annotation);
  cross-layer spearman(protein logFC, dTE).
- results: L1g_tests.tsv. Verdict after user review; honest reading: protein layer at 6h still integrates
  slow abundance dynamics; stage-delay confound documented by Li themselves (56/78 down = merozoite proteins).
- claim impact: direct protein-layer evidence for/against Asn axis under starvation; no claim change yet.


## Verdict (2026-09-04): PARTIAL POSITIVE at protein layer with the correctly-designed dataset
- spearman(Asn content, protein logFC under 6h starvation) = +0.167 (p=1.4e-8, n=1146);
  survives stage-delay control (excl. schizont-peak genes: +0.134, p=1.7e-5, n=1019).
- Distributional tilt, NOT hit-driven: decile MW p=0.11-0.17 n.s.; sig-set fractions at baseline.
- Cross-layer gene-by-gene (protein vs dTE): n.s. (expected - 6h abundance is a slow integral of TE).
- Contrast with D008 TMT (D-029, n.s.): the difference is DESIGN (starvation 6h, same-study TE sibling)
  vs (DHA 12h, 1249 proteins). Protein-layer support for the Asn axis exists but is shallow + aggregate.
- researcher note: PXD056396 was catalogued as 'tRip-KO' in the subagent brief; it is in fact the SAME
  iScience study as GSE226632, whose Table S2 IS the AA-starvation proteomics. Correction logged (D-031).
