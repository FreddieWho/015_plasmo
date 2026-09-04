# Session Handoff — G0 → M1

- session_id: 2026-09-03-G0-close
- date/time: 2026-09-03 UTC
- agent/model: Kimi Code CLI
- current milestone/gate: M1 / GATE_A (G0_READINESS → GO)
- objective: Close G0 readiness gate and start WP1 composition map (M1-01..04), no offset from 02/03/04 contract; honor 05/06 download gates
- files read: docs/M0_READINESS_AND_DATA_AUDIT.md v1.0; docs/manifests/project_state.yaml, gate_status.tsv, data_registry.tsv, frozen_assembly_candidates.tsv; data/checksums/* (6 manifests); data/veupathdb-file-download.zip (147M sha8683515); data/raw/ncbi-datasets/genomes/* (18 zips 454M), PlasmoDB-71 46 files
- tasks completed: Verified A-level D001-D005 DOWNLOADED_VERIFIED (veupathdb re-download OK, NCBI 18 genomes OK, taxonomy 18/18, UniProt 16 fasta, NAR 284 ERS); disk 2.7G <20G budget; flipped gate_status G0 REVIEW_REQUESTED→GO and project_state M0→M1 (GATE_A)
- outputs and paths: docs/manifests/gate_status.tsv (GO 2026-09-03); docs/manifests/project_state.yaml (M1_IN_PROGRESS, WP1, GATE_A); this handoff
- QC status: checksums OK (veupathdb, ncbi-genomes, taxonomy, uniprot, nar2017); species-audit SHA path-caveat noted but files verified via DOWNLOAD_MANIFEST.json
- failed attempts and why: species-audit sha256sum -c fails when run from repo root due to relative path mismatch (files are in data/raw/ncbi-datasets/species-audit-20260903/, not CWD); not a data integrity failure
- data registry changes: G0 evidence_refs appended with user re-download confirmation; D001-D005 remain DOWNLOADED_VERIFIED 2026-09-03
- evidence IDs added/updated: none new biological claims (M0 is data-audit only per 04 §10)
- claim IDs affected: none (DUAL_OBSERVATION_ONLY_UNTIL_GATE_B)
- gate implication: G0 GO unlocks WP1 M1-01..04; GATE_A still requires ≥2 independent composition contrasts with LCR-masked/unmasked dual-track, phylogenetically controlled, leave-one-clade-out
- hard blockers: None for M1; B-level (D006-D025) remains NOT_STARTED by design — no auto-download until GATE_A/B per 05 §4-§6
- user/expert decision needed: (1) P. vivax P01 vs Salvador-I sensitivity plan, (2) ovale downgraded presentation, (3) second deep-dive lineage priority (Pv/Pk/Pb) for expert input at Gate A review
- exact next task: M1-01 composition map — run data/derived/WP1/M1-01_composition/run_composition.py on frozen 18 (focus 12 FROZEN + 2 downgraded annotation-aware) to produce genome/compartment/CDS/GC1-3/4D/dinucleotide/homopolymer tables
- next task acceptance criterion: M1-01 outputs include per-species GC, N50, N-fraction, CDS GC/GC1-3/4D, dinucleotide, homopolymer, params.yaml + input_manifest.tsv + README per 04 §10; LCR dual-track noted; Pf mitochondrion/apicoplast absence flagged as method limit, not silent imputation
