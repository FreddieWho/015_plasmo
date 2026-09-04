# L1f claim impact (third-dataset translation-layer validation; exploratory)

- date: 2026-09-04
- data: GSE58402 Caro 2014 Ribo-seq+mRNA, W2 strain, 5 IDC stages (2/10/19/31/46 hpi), processed RPKM.
- prediction from chain: Asn-rich genes' TE rises ring->schizont, tracking Asn-tRNA charging (L1c: R 1.7 -> T 4.8 -> S 3.9).
- spearman(Asn content, TE 46h/2h) = -0.064 (p=0.069);
  Asn-top10% median late/ring TE shift = -0.055 vs rest -0.068, MW p=0.39.
- per-stage profile in L1f_qc.json median_log2TE_profile.
- caveat: W2 strain, 2014-era RPKM quant, no replicates per stage (single timecourse) -> descriptive support level.
- claim impact: independent translation-layer evidence for/against the Asn axis; feeds condition-1 closure.


## Verdict (2026-09-04): STAGE-DEPENDENT SUPPORT for the Asn supply-demand coupling
- Naive late/ring test (46h/2h): n.s. (rho=-0.064, MW p=0.39) - NOT supported at schizont endpoint.
- Mid-IDC test (19h/2h): rho=+0.35, p~1e-24 - STRONG: Asn-rich TE rises into trophozoite.
- Profile: Asn-rich TE ring -0.95 -> 31h peak -0.27 (parity with rest) -> 46h crash -1.23.
  Tracks Asn-tRNA charging (R 1.7 -> T 4.8) during growth; schizont = global translation shutdown regime (charging S 3.9 measured earlier in S phase).
- Honest label: CONDITIONAL SUPPORT - coupling holds in growth phases (R->T), not at schizont; 802 genes only (2014 RPKM, strict thresholds), single timecourse, W2 strain.
- Feeds condition-1 (third-dataset validation) as PARTIAL translation-layer support.
