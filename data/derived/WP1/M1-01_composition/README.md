# M1-01 Composition Map

**Scope** WP1/M1 Gate A descriptive composition only (04 §3 A-level). No translation axis yet (active_axis DUAL_OBSERVATION).

**Input** 18 frozen_accessions from `data/metadata/frozen_assembly_candidates.tsv` (T01). NCBI Datasets unpacked `genomic.fna` + `cds_from_genomic.fna` (D001). All A-level (05).

**Method**
- Genome GC/N/N50 via str.count (ATGC) — avoids per-base loops.
- Compartment via header.lower(): apicoplast/mitochondrion first; `chromosome Unknown` → nuclear_other; `chromosome` → nuclear_chromosome; `contig/scaffold/region` → nuclear_contig; else nuclear_other.
- CDS GC/GC1/GC2/GC3 via `seq[0::3]` slicing counts; GC4D = third base of 8 fourfold codon families (GCN/CGN/GGN/CTN/CCN/TCN/ACN/GTN) counting only third base.
- Dinucleotide 16 frequencies via numpy vectorized coding (A0 C1 G2 T3, N-masked) overlapping per contig.
- Homopolymer ≥5 via regex `A{5,}|T{5,}|G{5,}|C{5,}` (C-level).

**Outputs**
- `M1-01_composition_table.tsv` — 18 rows, 51 cols (genome/CDS/GC1-3/4D/homopolymer + 5 compartment len/count/gc/n)
- `M1-01_compartment_breakdown.tsv` — per-species per-compartment
- `M1-01_dinucleotide.tsv` — 16 dinucs per species
- `M1-01_homopolymer.tsv` — A/T/G/C runs ≥5 per species
- `M1-01_summary.md` — sorted by genome GC

**Key observations (descriptive, no claim)**
- Genome GC 18.2% (P. gaboni) → 52.3% (T. gondii), continuous gradient.
- SP009/SP010 (P. ovale) unannotated → CDS NA (expected, not error).
- Mito missing in Pf 3D7 / reichenowi / vivax etc. and api missing in some Pf relatives = NCBI assembly did not include organellar contigs (method limitation; do not fetch until Gate B per 04 §3).
- Homopolymer density anti-correlates with GC (Laverania high).

**Artifacts**
- `params.yaml`, `input_manifest.tsv` (sha256 per input/output), `run_composition.py` (302 lines).

**Reproduce**: `python3 data/derived/WP1/M1-01_composition/run_composition.py`
