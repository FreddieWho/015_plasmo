# M3-02 functional layer (D008 extraction + D009/D010 processed)

Q: does an independent functional layer (tRNA reprogramming + K13 perturbation + DHA transcript response) support the M2 tiered candidates (H3/CLM04 supporting only)?

## Inputs

- data/raw/pmc/PMC11153160/PMC11153160.pdf (31pp) + supp (29pp), DOI 10.1038/s41564-024-01664-3
- data/raw/geo/GSE151189/GSE151189_series_matrix.txt.gz (8.4MB, GPL18893, 156 cols: 5 K13 lines × time × ±DHA × 3rep)
- PRIDE PXD043747: metadata only (39 files all RAW/PEAK, no processed quants) — see PRIDE_PXD043747_NOTE.md

## Outputs

- M3-02_tRNA_mods.tsv: reprogrammed mods with direction + test + pdf ref
- M3-02_codon_bias.tsv: LysAAA/AAG driver stats + K13 52/57 AAA
- M3-02_assays.tsv: DHA/RSA/cKD/drug-panel/heat-shock conditions + readouts
- PRIDE_PXD043747_NOTE.md: no-raw-MS decision record

## Boundaries

- PDF text via pypdf; modification/codon strings verified against extraction (s2U/U34 naming kept as printed).
- Transcript abundance is NOT translation rate; TMT layer is paper-reported, not re-quantified here.
- No wet-lab proposals, no target ranking (Gate D).
