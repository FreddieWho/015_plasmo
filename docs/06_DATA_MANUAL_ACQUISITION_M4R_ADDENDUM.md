# M4R manual-acquisition addendum (NEW items only; does not repeat 06 §M01–M06)

**User rules:** official pages only; keep original filenames; place under `data/manual_inbox/<item_id>/`;
record date + URL + version; `sha256sum <file> > <file>.sha256`; notify agent for verification.
Agent then moves verified copies to read-only raw + registry (no author requests, no reviewer creds).

## T3 — manual / click-through / accession-unverified (do NOT auto-fetch by guessed number)

- **M4R-T3-01 — 2026 isogenic K13 scRNA DHA-6h (C-acute, highest gain): RESOLVED 2026-09-05.**
Front. Cell. Infect. Microbiol. 16 (2026), doi `10.3389/fcimb.2026.1816105` (Brief Research Report,
K13 C580 WT vs K13 580Y isogenic, 6 h DHA pulse, Seq-Well scRNA). Data availability verified from
article text: SRA BioProject **`PRJNA1049964`** (raw reads + processed UMI count matrices) +
Zenodo **`https://zenodo.org/records/20344254`** (matrices + code). Reference Pf 3D7 PlasmoDB v68.
Download: Zenodo record first (processed matrices, small), SRA only if re-alignment needed.
Target: `data/manual_inbox/M4R_K13_SCRNA_2026/` (agent: verify on arrival, then move to read-only
raw + registry). Note: transcript-only — does NOT close the translation-layer gap (M4R-X01 stands).
- **M4R-T3-02 — Pf8 release bundle (if automated fetch hits signed-URL/terms click):**
from https://www.malariagen.net/resource/36/ → processed VCF/Zarr/metadata/CNV subset only
(Pf, release version recorded; per-file list at download time). No FASTQ.
Target: `data/manual_inbox/M06_MALARIAGEN/` (existing convention; new subfolder per release).
- **M4R-T3-03 — GCN5 2026 supplementary (if publisher blocks direct fetch):**
from PMID 42323337 landing page → repeat-deletion + Py-complement processed tables (exact xlsx names
recorded). Target: `data/manual_inbox/M4R_GCN5_2026/`.
- **Per-species file note:** all T3 items above are single-species (Pf unless stated); multi-species
expansion (Pv/Pk/Pb ortholog checks) only on main-agent instruction with per-species file lists.

- **M4R-T3-04 — essentiality processed tables (D-boundary, T2→manual 2026-09-05):** automated fetch
BLOCKED — PMC serves supplements behind a JS proof-of-work interstitial (`cloudpmc-viewer-pow`, curl
gets 1.8 KB HTML); science.org is Cloudflare-walled; PkEssenDB is a Dash app (no static export URL).
User browser click-through (~5 min), keep original filenames, target
`data/manual_inbox/M4R_ESSENTIALITY/<paper>/` + sha256 + notify agent:
  - Elsworth 2025 (Science eadq6241, PMC12104972): `https://pmc.ncbi.nlm.nih.gov/articles/PMC12104972/`
    → Data_S1–S10.xlsx + Supp_Figures.pdf (`/articles/instance/12104972/bin/NIHMS2082619-supplement-*`).
    Alternative: PkEssenDB Download button `https://umbibio.math.umb.edu/PkEssenDB/`.
  - Oberstaller 2025 (Science eadq7347, PMC12131478): `https://pmc.ncbi.nlm.nih.gov/articles/PMC12131478/`
    → table_S1–S11.xlsx + 12.pdf (`/articles/instance/12131478/bin/NIHMS2081632-supplement-*`).
  - Zhang 2018 (Science aap7847, PMC6360947): `https://pmc.ncbi.nlm.nih.gov/articles/PMC6360947/`
    → Table_S1–S9.xlsx (`/articles/instance/6360947/bin/NIHMS1004827-supplement-*`).

## TX — do not pursue (explicit)
Author-request / lead-contact-only / private Drive-Dropbox / reviewer credentials / EGA-dbGaP
human-controlled / full FASTQ-WGS bulk / whole-AlphaFold-DB bulk / guessed-accession downloads.
(Zenodo public records such as 20344254 are allowed and encouraged.)
