# M4R-D ROUND2 REPORT (subagent, bounded per M4R-D_ROUND2_PACKET.md)

**Verdict: CONDITIONAL-STRONG for D-backbone.** Figure-3 verdict: **SURVIVES (frame-internal)**.
D2b keeps the architecture claim alive within the quantified sampling frame; D6 formal meets the
≥2-independent-transition rule (T2+T3); branch-LOO holds in all 7 (OR 1.71–4.35); D3 upgraded
where matchable; Pf8 boundary supports population-layer separation. Association only — no
mechanism claimed. T1 negative kept as registered negative, not rescued.

## Inputs (all in hand, no downloads; versions)
- Round-1 M4R-D outputs + scripts (run_m4rd.py D1/D2/D3/D5 functions replicated; round3 reuse
  pattern for xp2pf/pid2og/PROT/mafft). L2v2 site table (29,558 sites), 175 pass OGs,
  ML supermatrix tree, M1 GC states, grantham refs, NCBI GFF GCF_000002765.6 (xp2pf n=5354,
  IPR genes n=3622; =Pf3D7 by PF3D7_ locus_tag construction).
- kmerRBH wide orthogroups.tsv (D6 candidate enumeration + OG-age classes).
- Pf8 Zenodo 7 files (~14MB): 24,409 samples; 34 marker columns; FWS; CNV calls.
- Zhang Table_S5 (header row=1; per-gene MIS parsed), Elsworth/Oberstaller unused beyond round-1.
- mafft ~/.conda/envs/sc/bin/mafft --auto (LD_PRELOAD unset per round-1 pattern).
- New scripts: run_m4rd_r2.py (LOO/D2B/D3M/PF8 steps), run_m4rd_d6.py (D6 formal). Seed 20260905.

## D2b whole-proteome attribution — SURVIVES (frame-internal)
- Coupled-harboring Pf genes: 91 (from 91 coupled OGs via OG→SP001→PF3D7 mapping).
- Harbor vs length-decile×OG-age matched background (8 strata): LCR frac med 0.0 vs 0.0,
  MW p=0.10, strat-diff 0.0 → NOT more LCR-rich. Asn frac 0.064 vs 0.119, p=3.9e-33
  (strat-diff −0.031); polyN_max 2 vs 3, p=2.4e-18; n_IPR 5 vs 1, p=5.7e-28.
  → coupled turnover lives in structured, domain-rich, Asn-depleted conserved proteins.
- Frame bias quantified (core-181 members vs rest of proteome): length 259 vs 472 aa
  (p=7.8e-23); Asn 0.056 vs 0.120 (p=4.9e-75); LCR med 0.0 vs 0.0 (p=0.022).
  → "churn not repeat-confined" is frame-internal; stated as such. DOWNGRADED not triggered
  because harbor proteins show zero LCR excess under matching.
- LCR caller dual-track inherited from round-1 (convention + 2 sensitivities in D2 tables).

## D6 formal natural experiment — 2/3 transitions composition-consistent (rule MET)
- Pre-registered: T1_AT_Lav A={SP001,SP002} vs B={SP003} (expect A>B);
  T2_GC_vivax A={SP004-007} vs B={SP008,SP011-013} (expect A<B);
  T3_rodent A={SP011-013} vs B={SP008} (polarity-checked at runtime: 0.225 vs 0.244 → expect A>B).
- Coverage: 39 wide OGs aligned (AP2 4 with nsp≥8; PUF_RNA 27; CHROM 6 nulls); 106/117 rows OK,
  11 UNRESOLVED (mafft/coverage).
- T2: AP2 4/4 expected sign (p=0.0625), PUF_RNA 21/27 (p=0.0030), median diff −0.017.
  T3: AP2 3/4, PUF_RNA 18/25 (p=0.022), median diff +0.008/+0.018.
  T1: AP2 2/4, PUF_RNA 9/25 — NEGATIVE as registered (transition polarity tiny: 0.19 vs 0.18 GC;
  weak transition by construction, disclosed).
- CHROM null flat on all 3 transitions. Rule (≥2 independent) met by T2+T3.

## Branch-LOO (D1 pipeline faithful rerun ×7, pruned tree + recomputed ancestral GC)
- Coupled-vs-uncoupled Fisher OR (full run: 2.25, p=9.5e-10):
  laverania-drop 2.16 (p=1.9e-07); vivax-drop 1.71 (p=0.0065); malariae-drop 2.39;
  rodent-drop 2.41; avian-drop 2.74; piroplasm-drop 4.35; coccidian-drop 2.86.
  Significant in 7/7. Weakest under vivax-drop (high-GC pole removed — expected).
  Piroplasm-drop informative threshold jumps to 0.142 (dGC compression after removing
  high-GC outgroups) — disclosed, does not break significance.

## D3 multidim — partial upgrade, scope renamed
- Domain-shared background where matchable: K13/Kelch 531 sites @4.5% coupled vs known 0%;
  MDR1/ABC 5056 sites @2.0% vs known 0% (broad domain → weak restriction, disclosed);
  CRT/DHFR/DHPS: zero domain-shared core sites (families absent from core set).
- Within-OG control structurally impossible: all 14 known-gene OGs absent from core
  alignments (fast-evolving drug genes fail QC) — consistent with D6 frame finding, recorded
  as limitation not failure.
- Zhang MIS parsed for all 5 genes (0.12–0.13, dispensable-bin; no MIS contrast among them →
  matched-essentiality analysis uninformative, deferred honestly).
- Scope: "conservation-matched + domain-shared where available" (was "matched").

## Pf8 population boundary — supports D3 separation at population layer
- 34/34 resistance markers segregating (call rate ~1.0; strict single-clone FWS≥0.95 subset
  consistent): crt76 derived 0.56/0.64, dhfr51 0.83/0.81, dhfr108 0.95/0.95, dhfr164 0.21/0.26.
  → adaptation segregates in present-day populations at cross-species-conserved sites.
- CNV locus rates: GCH1 0.28, PM2/PM3 0.28, MDR1 0.20, CRT 0.086; HRP2/3 deletions ~0.45.
- SNP-level burden DEFERRED (Zarr streaming beyond round-2 budget) — stated, not silently dropped.

## QC / failures
- LOO_laverania first crashed (KeyError -1: stale parent pointer after tree-collapse);
  fixed by repairing parent on collapse; rerun passed. Lesson logged.
- D6: 2 proteomes missing from manifest (SP009/SP010 ovale GCA files) — skipped with log;
  all 16 analysis species loaded.
- D2B first run: IPR parsed from gene lines only (=0 genes); fixed to CDS+gene lines (3622).
- Zhang S5 header row=1 (not row 0); fixed after inspection.
- MWU always asymptotic; degenerate groups skipped+logged (D6 UNRESOLVED rows).
- No downloads; no AlphaFold/DL; M4R-X01 untouched; terminal layer kept.

## Claim impact
- Supports CLM08 (partition + axis-separation now with LOO-robust branch layer, frame-quantified
  architecture, population-layer D3 support). Weakens nothing; T1 negative bounds D6.
- Supports CLM09 D5-component (round-1 enrichment untouched).
- Prohibited extrapolations respected: no mechanism, no dual-axes, no repeat-first,
  D2b kept frame-internal, T1 not rescued.
- Suggested main-agent ledger: EVID-M4R-D-002 (ROUND2); update M5 convergence L4→full-ish
  pending A-round2 (L4 fully needs A-state too).

## Output paths (all under /home/huyudi/015_plasmo/data/derived/WP4R/M4R-D/)
M4RD_D1_LOO.tsv; M4RD_D2b_wholeproteome.tsv; M4RD_D2b_matched.tsv; M4RD_D2b_framebias.tsv;
M4RD_D2b_verdict.txt; M4RD_D3_multidim.tsv; M4RD_D4_Pf8boundary.tsv; M4RD_D4_Pf8cnv.tsv;
M4RD_D6formal_pertransition.tsv; M4RD_D6formal_summary.tsv; run_m4rd_r2.py; run_m4rd_d6.py;
claim_impact.md (appended); input_manifest.tsv (appended); checksums.sha256 (regenerated).
