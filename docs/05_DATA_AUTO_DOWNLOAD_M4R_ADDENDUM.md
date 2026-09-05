# M4R auto-download addendum (NEW items only; does not repeat 05 §A–D)

**Rule:** processed-first; ≤20 GB auto / 20–100 GB Gate-gated / >100 GB stop-and-ask. Re-verify every
accession at fetch; record DOWNLOAD_MANIFEST.json + sha256; update data_registry (no Gate/claim edits here).

## T1 — fetch processed matrix/metadata now (small, direct)

- **M4R-A02 GSE75795:** GEO FTP `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75795/`
→ series matrix + SOFT + metadata. Small. No SRA unless Gate-justified.
- **M4R-A03 ChIP trio GSE120448 / GSE134268 / GSE120488:** per-accession GEO FTP
(`.../GSE120nnn/<ACC>/`, `.../GSE134nnn/GSE134268/`) → processed peaks/matrices + metadata.
- **M4R-C01 GSE225340:** GEO FTP `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE225nnn/GSE225340/`
→ series matrix + metadata. Small-to-moderate.
- **M4R-C02b GSE59099:** GEO FTP `.../GSE59nnn/GSE59099/` → series matrix + metadata (chronic
transcript end; Mok protein end already in-house).
- **M4R-D02 essentiality:** PkEssenDB export (https://umbibio.math.umb.edu/PkEssenDB/) + PMC
supplements for Elsworth/Oberstaller/Zhang → processed tables only. Small.

## T2 — conditional (high value, resolve names at fetch, no bulk raw)

- **M4R-D01 Pf8:** from https://www.malariagen.net/resource/36/ → processed VCF/Zarr/metadata/CNV
subset only (record release version + exact filenames at fetch). Never FASTQ.
- **M4R-A01 GCN5 2026:** publisher/PMC supplementary xlsx for PMID 42323337 (full DOI re-verify;
record exact filenames). If landing page exposes direct xlsx links, treat as T1.
- **PXD043916 tRip-KO:** PRIDE `https://www.ebi.ac.uk/pride/archive/projects/PXD043916` → verify
organism/files first; processed results only; halt if organism ≠ Plasmodium or files mismatch.
