# Manuscript supplementary tables — provenance index (A/C-side lane)

All tables are verbatim copies of frozen analysis outputs (no recomputation).
D-side tables (TableS1 D1 concordance, TableS2 D2b, TableS3 D6) are delivered by
the FIG-D lane; see its report on arrival. Do not duplicate them here.

| Table file | Source (frozen) | Source sha256 | Rows | Notes |
|---|---|---|---|---|
| TableS4_interpro_enrichment.tsv | data/derived/WP4R/M4R-A/M4RA_interpro_enrichment.tsv | see data/derived/WP4R/M4R-A/checksums.sha256 + checksums_round2.sha256 | 17 + header | formal-vs-regex side-by-side; sexual formal == AP2 set (no double-count); PUF strict n=2 |
| TableS5_lifecycle_tripole.tsv | data/derived/WP4R/M4R-A/M4RA_lifecycle_tripole.tsv | as above | 28 + header | gametocyte/zygote/liver; liver EXCLUDED per detection aux |
| TableS6_C_partition_asnfrac.tsv | data/derived/WP4R/M4R-C/M4RC_acute_chronic_partition.tsv | see data/derived/WP4R/M4R-C/checksums.sha256 | 142 + header | asn_frac rows; q column merged from M4RC_validation_table.tsv in Fig5a panel table |

Figure panel tables live in manuscript/figures/panels/ (fig4a_forest, fig4b_strata,
fig4c_gcn5_curves, fig4c_pbmap, fig4c_tripole_box, fig5a_partition) and are strict
subsets/derivations (Fisher-z CI only) of the frozen tables above.

| TableS1_D1_concordance.tsv | M4RD_D1_concordance_summary.tsv | M4R-D checksums.sha256 | 7 + header | coupled 0.791 vs uncoupled 0.628, OR=2.25 |
| TableS1b_D1_LOO.tsv | M4RD_D1_LOO.tsv | as above | 7 + header | branch-LOO OR 1.71–4.35, 7/7 significant |
| TableS2_D2b_matched.tsv | M4RD_D2b_matched.tsv | as above | harbor n=91 vs matched bg |
| TableS2b_D2b_framebias.tsv | M4RD_D2b_framebias.tsv | as above | core-181 vs rest sampling-frame quant |
| TableS3_D6formal_summary.tsv | M4RD_D6formal_summary.tsv | as above | T2+T3 rule met, T1 negative, CHROM null |
| TableS3b_D6formal_pertransition.tsv | M4RD_D6formal_pertransition.tsv | as above | per-transition sign tests |
| TableS3c_D3_multidim.tsv | M4RD_D3_multidim.tsv | as above | domain-shared where available |
| TableS3d_Pf8boundary.tsv | M4RD_D4_Pf8boundary.tsv | as above | 34/34 markers segregating + CNV |
