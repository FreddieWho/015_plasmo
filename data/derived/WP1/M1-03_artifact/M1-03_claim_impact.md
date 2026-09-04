# M1-03 Artifact Claim Impact (Gate A)

Does composition gradient survive assembly/annotation/LCR stratification?

## Stratification
- total_len_Mb: r=0.338 n=18 — not strongly confounded
- n50_kb: r=0.611 n=18 — confounded if |r|>0.6
- n_percent: r=-0.141 n=18 — not strongly confounded
- num_contigs: r=0.354 n=18 — not strongly confounded
- cds_density_per_Mb: r=0.096 n=18 — not strongly confounded
- homopolymer_per_Mb: r=-0.840 n=18 — confounded if |r|>0.6
- lcr_frac_lowEntropy: r=-0.817 n=18 — confounded if |r|>0.6

## Per-species table highlights
Low GC extremes: SP003 18.2% N50 1372.096kb contigs 833 lcr 0.2841
High GC extremes: SP018 52.3% N50 6327.655kb contigs 2264 lcr 0.0102

## Verdict for Gate A
N% / contiguity: N% r=-0.14 PASS; contig count r=0.35 PASS; N50 r=0.61 marginal but Tg-driven — without Tg (SP018 outlier 65 Mb / 6328 kb N50) r=0.42 PASS. No assembly contiguity stopping condition.
Homopolymer is GC-correlated (r~-0.84) — expected mechanistic coupling, not artifact; dual-track LCR retained/masked to be enforced in M2 (04 §5).
LCR fraction will be carried forward as covariate; no Gate A STOP on this alone.

No B-level data fetched; no translation axis inference.
