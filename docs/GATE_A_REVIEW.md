# Gate A Review — M1 Composition Evolution

- gate_id: GATE_A
- review date: 2026-09-03
- decision candidates: GO / PIVOT / STOP-DOWNGRADE / HOLD
- pre-specified criterion: Gate A — at least two relatively independent composition contrasts, compartment/polarity repeatable, not driven solely by annotation/LCR/assembly (03 §3)
- supporting evidence IDs: EVID-M1-001 (composition map), EVID-M1-002 (k-mer phylogeny + ancestral + independence), EVID-M1-003 (artifact screen)
- contradicting evidence IDs: none blocking; N50 marginal r=0.61 flagged but Tg-driven (without Tg r=0.42)
- data quality limitations: SP009/SP010 ovale kurtisi/wallikeri have no cds_from_genomic.fna (CDS GC NA, flagged NEED_REANNOTATION); 4 species lack api/mito contigs in NCBI assembly (SP001/003/012/013 api=0, method limitation not data error); 6 scaffold-level assemblies retained but FROZEN
- sensitivity results: leave-one-clade-out 7 clades all retain gradient (range 0.23–0.34, full 0.34); N%/contig/Len not strongly confounded; LCR/homopolymer strongly GC-correlated as expected mechanistic coupling (see M1-03)
- independent lineages/replicates: Grade A — independent low GC in Laverania (SP001-003 18–19%) and avian P. relictum (SP014 18.3%, non-sister); independent high GC in vivax/knowlesi clade (SP004-007 38–40%) and piroplasm (Babesia/Theileria 34–41%) and coccidian outgroup (Tg 52% separate) => ≥2 independent contrasts satisfied
- unresolved alternatives: NAR mutation-spectrum vs substitution distinction deferred to WP2 with explicit nulls; compartment-specific GC1-3/4D gradient qualitatively follows genome GC but formal AA-fixed codon analysis is Gate B; vivax P01 vs Salvador-I dual-run deferred per open_decisions
- minimum next information gain task (only for HOLD): n/a
- cost of next task: WP2 codon/AA decomposition uses existing CDS + orthogroups + LCR dual-track, no new raw data fetch
- expert input: requested for deep-dive second species (Pv vs Pk vs Pb) before WP2 codon design
- final decision: **GO**
- active axis after decision: **DUAL_OBSERVATION_ONLY_UNTIL_GATE_B** (unchanged; Gate B will select TRANSLATION_PRIMARY / PLASTICITY_PRIMARY / DOWNGRADE per 03 §4 hard rule)
- claims upgraded/downgraded: CLM01 HYPOTHESIS → SUPPORTED_CANDIDATE (composition transitions observed, independent, not artifact-only); CLM02-CLM07 remain HYPOTHESIS/FALLBACK (no codon/functional inference at Gate A)
- project_state updated: yes — M1_COMPLETE, Gate A GO 2026-09-03, WP1 evidence ledger closed

## Evidence Summary (WP1)

### M1-01 Composition Map (18 species, frozen accessions)
- Genome GC spans 18.2% (P. gaboni SP003) → 52.3% (T. gondii SP018), range 34.1%. Plasmodium-internal 18.2% → 41.6% (B. bovis). Ordering conserved across compartments: CDS GC tracks genome GC; GC3 and 4D GC show steepest slope (e.g. SP003 GC3 15.6%/GC4D 15.0% vs SP004 GC3 53.0%/GC4D 58.1%).
- Table: data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv; breakdown: M1-01_compartment_breakdown.tsv; dinucleotide + homopolymer + summary included.
- No B-level data used.

### M1-02 Phylogeny + Ancestral + Independence
- Method: k=21 stride-20 crc32 min-hash sketch 10k, Jaccard distance, Bio.Phylo NJ rooted on Tg SP018. Composition-independent (not GC).
- Tree: data/derived/WP1/M1-02_phylogeny/M1-02_species_tree.nwk (ascii in M1-02_tree_ascii.txt); distance in M1-02_kmer_distance.tsv.
- Ancestral GC (squared-change parsimony, mean children): 17 internal nodes in M1-02_ancestral_gc.tsv; root ~43.6% (includes Tg), Plasmodium-ancestor ~34.8%, Laverania+avian ancestor ~18.5%.
- Leave-one-clade-out (7 clades): all retain gradient (Laverania-removed still 18.3–52.3% range; coccidian-removed still 18.2–41.6%).
- Independence grade A: low GC in Laverania + avian (separate branches); high GC in vivax clade + piroplasm (separate).
- Note: distance values compressed near 1.0 for k=21 at this divergence — expected, tree is illustrative for independence claim only; full orthogroup tree can be built in WP2 if needed.

### M1-03 Artifact Screen
- Metrics vs genome GC (Pearson, n=18): total_len r=0.34, n_percent r=-0.14, num_contigs r=0.35, cds_density r=0.10 — not strongly confounded; n50 r=0.61 marginal but driven by Tg 65 Mb/6328 kb outlier (without Tg r=0.42).
- Homopolymer r=-0.84 and LCR low-entropy fraction r=-0.82 — strong negative correlation, interpreted as mechanistic coupling (AT-rich genomes generate poly-A/T and low entropy) not primary artifact; dual-track LCR retained/masked enforced as Gate B covariate per 04 §5.
- No assembly/annotation stopping condition. Details: data/derived/WP1/M1-03_artifact/M1-03_artifact_table.tsv + correlations + claim_impact.

### Gate A Verdict
GO — at least two relatively independent composition transitions are observed, span compartments, survive artifact stratification. WP1 satisfies 03 §3 and 04 §10 (DUAL_OBSERVATION). No claim beyond C1 is made.

## What Was NOT Done (per 05/06 gating)
- No D006–D025 fetched at Gate A (all B/C-level: reads, structures, host omics, drug evolution). D001–D005 already DOWNLOADED_VERIFIED satisfies M0–M1 A-level. Pausing B/C downloads until Gate B is compliant, not a deviation.
- No codon/AA residual, translation, or plasticity inference (Gate B/C).
- No ML ranking or target nomination.

## Next
WP2 will: AA-fixed synonymous nulls, local composition controls, orthologous-site/ancestral/structure/LCR dual-track, phylogenetically aware tests, leave-one-clade robustness — then Gate B axis selection.

## Pointers
- M1-01: data/derived/WP1/M1-01_composition/
- M1-02: data/derived/WP1/M1-02_phylogeny/
- M1-03: data/derived/WP1/M1-03_artifact/
- Frozen panel: data/metadata/frozen_assembly_candidates.tsv
- Data registry: docs/manifests/data_registry.tsv (D001–D005 VERIFIED)
