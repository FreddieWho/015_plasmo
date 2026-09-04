# M0 Readiness and Data Audit — PLASMODIUM-C2F

**Project:** PLASMODIUM-C2F (疟原虫组成—功能, Composition → Codon → Function)  
**Evidence standard:** Nature Microbiology (constraint, not outcome commitment)  
**Milestone:** M0 / WP0 — 立项与数据审计  
**Gate:** G0_READINESS (Are the data, species contrasts and public access sufficient to start?)  
**Date:** 2026-09-03 (UTC)  
**Version:** v1.0 — first deliverable; no biological analysis beyond data audit  
**Authors:** Agent + 18-species audit (NCBI Datasets 18.36.0, VEuPathDB Release 71, UniProt 2026_03, NAR/PMC5389722)  
**Active axis:** DUAL_OBSERVATION_ONLY_UNTIL_GATE_B (per `project_state.yaml`)  

> 约束：本报告只回答“哪些问题可以被哪些数据回答”，不做组成差异、密码子偏性或机制方向的提前锁定。所有下游分析需通过 Gate A/B/C/D 的独立证据门。

---

## 0. 执行摘要 / Go-No-Go

| G0 条件 | 阈值 | 本次审计结果 | 判定 |
|---|---|---|---|
| 高质量 Plasmodium 候选 | ≥8 | **14** 个 Plasmodium 中 12 个 FROZEN 可直接用于核心比较，另 2 个 (SP009/SP010 ovale) 为 annotated-scaffold 降级可用，总计 12–14 | ✅ PASS |
| 独立组成对照潜力 | ≥2 | 至少 **3** 个相对独立对照候选：Laverania 内 3D7 高AT vs 近缘种 / vivax/knowlesi 中等GC vs Laverania 高AT / rodent/avian/piroplasm GC梯度与 Laverania 的跨支对比（需 M1 以 phylogenetically controlled composition 正式检验） | ✅ PASS (candidate-level) |
| 深挖谱系数据 | ≥2 (Pf + 第二谱系) | Pf 深挖层就绪 (D005/D006-D013 预留)；第二谱系 P. vivax (Pv4), P. knowlesi (scRNA + Tn-seq), P. berghei (lifecycle scRNA) 均为公开可获取，M1 前不拉取 raw reads | ✅ PASS |
| 公开可获取性 | 无作者请求依赖 | D001–D005 全部 `DOWNLOADED_VERIFIED`，仅公共库；其余 M3/M4 层标注为 CONDITIONAL/OPTIONAL，按 Gate 触发下载 | ✅ PASS |

**总体判定：M0 具备启动 WP1 的数据基础，可申请 G0 评审后进入 M1。** 不建议在 G0 关闭前启动任何大规模群体 reads / 全量结构预测 / 宿主组学。

---

## 1. 执行包完整性验证

| 项 | 预期 | 实测 | 状态 |
|---|---|---|---|
| `docs/manifests/project_state.yaml` | package v1.0, scientific_spec 2026.09, M0/WP0/G0 | 已读：current_milestone M0, active_gate G0_READINESS, frozen_invariants 11 条完整 | ✅ |
| `docs/manifests/species_panel_seed.tsv` | 候选物种种子表 | 18 物种已审计；对应 taxid 全部 NCBI Taxonomy 验证通过 | ✅ |
| `docs/00–12` 执行包文档 | 12 份 + manifests + templates + sources | 13 份 md 已通读；`sources/` 含 iScience / NAR / Proposal 三篇起始工作 | ✅ |
| 三篇起始工作版本 | Proposal CN docx + NAR gkw1259 + iScience seed | 已定位，DOI 10.1093/nar/gkw1259 对应 PMC5389722 已下载校验 | ✅ |
| 存储约定 | `data/raw/<source>/<dataset_id>/` + `DOWNLOAD_MANIFEST.json` + `data/checksums/*.sha256` | 本次全部遵守（见 §4） | ✅ |

变更控制：`templates/T01_CHANGE_REQUEST.md`；任何 assembly 冻结变更需走 change request。

---

## 2. 候选物种当前参考组装审计与冻结建议

### 2.1 审计方法

- 工具：`NCBI datasets CLI 18.36.0` (`/tmp/datasets`) + NCBI Datasets API (`summary genome taxon`)
- 查询：18 个 taxid 各自 taxon 级别检索，记录 `total assemblies / GCF / GCA`, 再对每物种 `best accession` 做 `summary genome accession --as-json-lines` 抽取 `assembly_level, total_length, gc_percent, N50, annotation_status/provider/release_date, gene_counts`
- 输出：`data/metadata/species_assembly_audit.tsv` (19 行含表头, 18 物种) + `data/metadata/frozen_assembly_candidates.tsv` (19 行) + 原始 JSON `data/raw/ncbi-datasets/species-audit-20260903/SP*_*.json` (18 个) + `DOWNLOAD_MANIFEST.json`

### 2.2 冻结总表（18 物种，2026-09-03 当前参考）

| species_id | 物种 | taxid | 冻结 accession | 库 | 水平 | 注释 | 长度 bp | GC% | N50 (bp) | 基因总数 / 蛋白编码 | PlasmoDB-71 对应 | Gate | 冻结动作 | 风险旗 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| SP001 | *P. falciparum* 3D7 | 5833 | **GCF_000002765.6** | RefSeq | Complete | annotated (2020-09-10) | 23,292,622 | 19.5 | 1,687,656 | 5,615 / 5,282 | Pfalciparum3D7 ✓ | G0/A | FROZEN | LOW |
| SP002 | *P. reichenowi* CDC | 5854 | **GCF_001601855.1** | RefSeq | Chromosome | annotated (2017-02-17) | 20,477,502 | 18.5 | 1,412,118 | 5,255 / 5,040 | PreichenowiCDC ✓ | G0/A | FROZEN | LOW |
| SP003 | *P. gaboni* G01 | 647221 | **GCF_001602025.1** | RefSeq | Chromosome | annotated (2016-10-21) | 20,356,150 | 18.0 | 1,372,096 | 5,429 / 5,196 | NOT_IN_PlasmoDB-71 | G0/A | FROZEN_GCF | LOW |
| SP004 | *P. vivax* **P01** | 5855 | **GCA_900093555.2** | GenBank | Chromosome | annotated | 29,040,213 | 42.5* | 1,678,596* | 6,519 PC† | PvivaxP01 ✓ | G0/A | FROZEN_P01 | SEED_CONFLICT_RESOLVED |
| SP005 | *P. knowlesi* H | 5850 | **GCF_000006355.2** | RefSeq | Chromosome | annotated (2023-04-03) | 24,359,384 | 38.5 | 2,162,603 | 5,464 / 5,324 | PknowlesiH ✓ | G0/A | FROZEN | LOW |
| SP006 | *P. cynomolgi* B | 5827 | **GCF_000321355.1** | RefSeq | Chromosome | annotated (2013-04-25) | 26,181,343 | 40.5 | 1,717,921 | 5,776 / 5,716 | NOT_IN_PlasmoDB-71 | G0/A | FROZEN_GCF_WITH_WARNING | ANNOTATION_GAP |
| SP007 | *P. coatneyi* | 208452 | **GCF_001680005.1** | RefSeq | Chromosome | annotated (2017-02-02) | 27,685,530 | 39.5 | 1,997,354 | 5,575 / 5,516 | NOT_IN | G0/A | FROZEN | LOW |
| SP008 | *P. malariae* UG01 | 5858 | **GCF_900090045.1** | RefSeq | Chromosome | annotated (2023-04-03) | 33,577,742 | 24.5 | 2,312,277 | 6,677 / 5,926 | PmalariaeUG01 ✓ | G0/A | FROZEN | LOW |
| SP009 | *P. ovale curtisi* GH01 | 864141 | **GCA_060229945.1** ⚠️ | GenBank | Chromosome | **unannotated** | 36,010,039 | 28.0 | 46,061 (scaffold annotated 8825 genes) | — / — | PovalecurtisiGH01 (scaffold) | G0/A | CANDIDATE_DOWNGRADED | NEED_REANNOTATION |
| SP010 | *P. ovale wallikeri* | 864142 | **GCA_060229965.1** ⚠️ | GenBank | Chromosome | **unannotated** | 35,654,824 | 29.0 | 173,728 (scaffold annotated 8421 genes) | — / — | NOT_IN (sister) | G0/A | CANDIDATE_DOWNGRADED | NEED_REANNOTATION |
| SP011 | *P. berghei* ANKA | 5821 | **GCF_900002375.2** | RefSeq | Chromosome | annotated (2023-04-03) | 18,740,246 | 22.0 | 1,640,193 | 5,216 / 4,933 | PbergheiANKA ✓ | G0/A | FROZEN | LOW |
| SP012 | *P. yoelii* 17X | 5861 | **GCF_900002385.2** | RefSeq | Complete | annotated (2020-05-23) | 23,043,114 | 21.5 | 2,046,250 | 6,233 / 6,037 | Pyoelii17X ✓ | G0/A | FROZEN | LOW |
| SP013 | *P. chabaudi* AS | 31271 | **GCF_900002335.3** | RefSeq | Chromosome | annotated (2023-04-12) | 18,938,342 | 23.5 | 1,632,254 | 5,389 / 5,199 | NOT_IN | G0/A | FROZEN | LOW |
| SP014 | *P. relictum* SGS1 | 85471 | **GCF_900005765.1** | RefSeq | Chromosome | annotated (2023-04-03) | 22,571,969 | 18.5 | 1,287,098 | 5,301 / 5,132 | NOT_IN | G0/A | FROZEN | LOW |
| SP015 | *Babesia bovis* T2Bo | 5865 | **GCF_000165395.2** | RefSeq | Chromosome | annotated (2022-11-03) | 8,195,476 | 41.5 | 1,797,577 | 3,988 / 3,959 | NOT_IN (PiroplasmaDB) | G0/A | FROZEN_GCF | NEWER_UNANNOTATED |
| SP016 | *Babesia microti* RI | 5868 | **GCF_000691945.2** | RefSeq | Chromosome | annotated (2023-04-04) | 6,395,281 | 36.5 | 1,766,409 | 3,685 / 3,601 | NOT_IN | G0/A | FROZEN | LOW |
| SP017 | *Theileria parva* Muguga | 5875 | **GCF_000165365.1** | RefSeq | Chromosome | annotated (2023-12-06) | 8,308,027 | 34.0 | 1,971,884 | 4,128 / 4,046 | NOT_IN | G0/A | FROZEN | LOW |
| SP018 | *Toxoplasma gondii* ME49 | 5811 | **GCF_000006565.2** | RefSeq | Chromosome | annotated (2023-04-03) | 65,633,124 | 52.5 | 4,973,582 | 8,925 / 8,318 | NOT_IN (ToxoDB) | G0/A | FROZEN | LOW |

\* SP004 P01 的长度与 N50 取自 GCA_900093555.2 的 NCBI 记录；GCF Salvador-I (GCF_000002415.2) 为 27.0 Mb / 1.16 Mb 作对照。  
† P01 蛋白编码基因数 6,519 来自 PlasmoDB PvP01_v2 与 literature_anchor 对照，显著高于 Salvador-I 的 5,392 — 这是冻结 P01 的关键依据。

**Alternative accessions（仅作灵敏度对照，不加入主冻结）：**
- SP004 Salvador-I `GCF_000002415.2` (5,392 PC, 27.0 Mb) — 已记录，必要时可作 P01 vs Salvador-I 双跑以检验亚端粒差异是否由 assembly 驱动。
- SP006 `GCA_058584565.1` (37.0 Mb, N50 3.1 Mb, **未注释**) — 仅在完成重注释后才考虑晋升。
- SP009/010 annotated scaffold `GCA_900088565.1` / `GCA_900088545.1` — 作为 ovale 的 interim annotated 替代，但碎片化严重 (N50 46k/133k)，不得用于 Gate A 的主要“独立对照”断言。

**质量等级与使用边界（按 04_ANALYSIS_CONTRACT §2）：**
- 12 个 FROZEN (LOW)：满足 `Complete/Chromosome + annotated + N50 >1 Mb + 基因数与近缘种同量级`，可直接用于发现层与 M1 组成地图。
- 2 个 FROZEN_GCF / FROZEN_GCF_WITH_WARNING (SP003/SP006/SP015)：可用但需在 manuscript 中声明 alternative 与 annotation gap；SP006 的较新未注释 assembly 不得静默替换。
- 2 个 CANDIDATE_DOWNGRADED (SP009/SP010)：染色体水平但无注释，scaffold 注释碎片化；仅作探索性对照，M1 正式 Gate A 需标注 `NEED_REANNOTATION`，若要晋升需运行可重复注释流程或等待 RefSeq 更新。

---

## 3. VEuPathDB Release 71 与 NCBI 冻结的交叉校验

### 3.1 PlasmoDB-71 落地

- **来源：** `data/veupathdb-file-download.zip` 147 MB，`sha256 86835150bc8e4166a64d1e967f089a18a280672b90a646ea4c047bfeb8953edf`
  - 备注：2026-09-02 初次下载仅 8 MB 且 `unzip` 报“中央目录缺失”（截断）；已按 05 重新下载校验通过。
- **解压：** `data/raw/veupathdb/PlasmoDB-71/` — **46 files, 554 MB**，涵盖 10 物种 × 5 文件类型（Genome.fasta / AnnotatedProteins.fasta / AnnotatedTranscripts.fasta / GFF / GeneAliases.txt，另含部分 InterPro GFF 子集）
  - 物种清单：`PbergheiANKA / Pfalciparum3D7 / PinuiSanAntonio1 / PknowlesiH / PmalariaeUG01 / PovalecurtisiGH01 / PpraefalciparumG01 / PreichenowiCDC / PvivaxP01 / Pyoelii17X`
- **清单与校验：** `data/raw/veupathdb/PlasmoDB-71/DOWNLOAD_MANIFEST.json`（46 条 per-file sha256，总计 553,765,852 bytes）+ `data/checksums/veupathdb-PlasmoDB-71.sha256` (`sha256sum -c` OK)

### 3.2 与 NCBI 冻结的一致性

| 维度 | PlasmoDB-71 | NCBI 冻结 | 一致性 |
|---|---|---|---|
| *P. falciparum* 3D7 | Pfalciparum3D7 (Pf3D7_v3) | GCF_000002765.6 (3D7) | ✅ 同株；基因数与 GC 一致 |
| *P. reichenowi* CDC | PreichenowiCDC | GCF_001601855.1 CDC | ✅ 同株 |
| *P. vivax* | **PvivaxP01** (PvP01_v2) | **GCA_900093555.2 P01** (冻结) vs GCF_000002415.2 Salvador-I (对照) | ✅ 种内株级分歧已显式记录：PlasmoDB 选用 P01，NCBI RefSeq 默认 Salvador-I；M0 选择与 PlasmoDB 一致的 P01，保留 Salvador-I 作灵敏度分支 |
| *P. knowlesi* H | PknowlesiH (PKNH_v2) | GCF_000006355.2 H | ✅ |
| *P. malariae* UG01 | PmalariaeUG01 | GCF_900090045.1 UG01 | ✅ |
| *P. ovale curtisi* GH01 | PovalecurtisiGH01 | GCA_060229945.1 (chr unannot) / GCA_900088565.1 (scaffold annot) | ⚠️ PlasmoDB GH01 基于旧 scaffold；NCBI 新 chr 未注释 — 不得混用 release，需在 M1 标注“跨库版本差” |
| *P. berghei* ANKA | PbergheiANKA | GCF_900002375.2 | ✅ |
| *P. yoelii* 17X | Pyoelii17X | GCF_900002385.2 | ✅ |
| *P. gaboni / P. cynomolgi / P. coatneyi / P. relictum / Babesia / Theileria / Toxo* | **不在** PlasmoDB-71（分别归属 PiroplasmaDB/ToxoDB 或新种） | 均已通过 NCBI 冻结 | ✅ 已记录“NOT_IN_PlasmoDB-71”，后续不强行对齐；必要时到对应 VEuPathDB 子库验证 |

**M1 使用规则：** 同一物种不得混用不同 release 的 genome 与 GFF；若同时使用 NCBI 与 PlasmoDB，需在 Figure 1 审计图中并排展示两者 GC/基因数/N50 差异，并对关键候选做双注释交叉检查（见 04 §2）。

---

## 4. 数据可访问性、落地路径与磁盘预算

### 4.1 A 级底座（D001–D005）— 已全部 `DOWNLOADED_VERIFIED`

| dataset_id | 来源 | 落地路径 | 内容与规模 | 校验 |
|---|---|---|---|---|
| **D001** | NCBI Datasets | `data/raw/ncbi-datasets/genomes/*.zip` (18 zips, **454 MB**) + `data/raw/ncbi-datasets/<ACC>/` (18 解压目录，合计 **1.5 GB** unpacked) + 审计 `data/raw/ncbi-datasets/species-audit-20260903/` (18 JSON) | 每 accession 含 `genomic.fna / cds.fna / protein.faa / genomic.gff / gbff / seq-report`；`species_assembly_audit.tsv` + `frozen_assembly_candidates.tsv` 各 19 行 | `data/checksums/ncbi-genomes-20260903.sha256` (18 条) `sha256sum -c` OK；`species-audit-20260903.sha256` OK；每目录 `DOWNLOAD_MANIFEST.json` 含 tool version / URL / size / sha256 |
| **D002** | VEuPathDB PlasmoDB-71 | `data/raw/veupathdb/PlasmoDB-71/` (46 files, 554 MB) | 10 物种 × Genome/Proteins/Transcripts/GFF/GeneAliases | `DOWNLOAD_MANIFEST.json` (46 条) + `veupathdb-PlasmoDB-71.sha256` OK |
| **D003** | NCBI Taxonomy + 文献生态 | `data/raw/ncbi_taxonomy/` (54 xml/json) + `data/metadata/taxonomy_host_ecology.tsv` (19 行) + `data/metadata/species_assembly_audit.tsv` | 18/18 taxid 验证通过：5833, 5854, 647221, 5855, 5850, 5827, 208452, 5858, 864141, 864142, 5821, 5861, 31271, 85471, 5865, 5868, 5875, 5811；生态位按 Laverania/vivax-clade/rodent/avian/piroplasm 分支标注 | `data/checksums/taxonomy-20260903.sha256` OK |
| **D004** | UniProt (host hemoglobin) | `data/raw/uniprot/hemoglobin/` (16 fasta) + `data/metadata/uniprot_hemoglobin_manifest.tsv` (17 行) | 8 宿主 × HBA/HBB (human, chimp, macaque, rodent, bird, cattle, broad host 组合)；UniProt 2026_03 release；mouse/chicken 各 2 entries，*M. fascicularis* HBA 为 unreviewed fallback | `data/checksums/uniprot-hemoglobin-20260903.sha256` OK；manifest 含 gene/query/resolved_url/release/size/entries/sha256 |
| **D005** | NAR 2017 / PMC5389722 / ENA | `data/raw/nar2017/` (70 MB) + `data/metadata/nar2017_ena_accessions.tsv` (285 行, 284 ENA) | PMC5389722 supplementary zip 36665988 bytes, sha `0af731ca…`；3 files 解压；Table S2 解析 284 ENA accessions ERS571523–ERS572179 (含 Bam/Clone/≥5/≥10 列) | `data/checksums/nar2017-20260903.sha256` OK；manifest `data/raw/nar2017/DOWNLOAD_MANIFEST.json` |

> **D001–D005 全部为 public_direct / public_api / public_repository，无需作者请求，符合 06 “No author-request-only”。**

### 4.2 存储审计

```
104K  data/raw/uniprot
296K  data/raw/ncbi_taxonomy
 70M  data/raw/nar2017
529M  data/raw/veupathdb
2.1G  data/raw/ncbi-datasets          # 包含 454M zips + ~1.5G unpacked + 800K audit JSON
----------------------------------------
2.7G  data/raw/ TOTAL  (2026-09-03 03:46 UTC)
 28K  data/checksums/ (6 manifests)
 68K  data/metadata/ (5 tsv)
```

- **预算：** 当前 <3 GB，远低于 docs 要求 **>20 GB 需用户确认** 的阈值，允许在 G0 后按需增量下载 D006–D013 处理后矩阵。
- **增量预估：** D006–D013 (GEO/PRIDE/MCA 处理后表格与伪 bulk) 预计额外 <5 GB；若未来触发 D016 raw reads (>100 GB) 必须单独审批并走 `T01_CHANGE_REQUEST.md`。

### 4.3 待触发层（D006–D025）— 审计可访问性但暂不下载

| 层 | 数据集 | URL | 可访问性预检 | M0 动作 |
|---|---|---|---|---|
| DEEP M3 翻译轴 | D006 GSE226632, D007 PXD056396, D008 PMC11153160, D009 PXD043747, D010 GSE151189, D011–D013 MCA | GEO/PRIDE/PMC/MCA 公开 | 均为公开 accession，无需新权限；PRIDE 需按“processed first”仅取 quantification/design | **不下载**（Gate C 前禁止全量 raw） |
| FUNCTIONAL M3/M4 | D015 Luth DHODH xlsx (PMC11809290), D017 Zenodo 13766513, D018 PkEssenDB PRJNA1116591 | PMC/Zenodo/UMB 公开 | 7 个 supplement xlsx + code archive 均为公共 | **不下载**，M3 触发时再取处理后表格 |
| POPULATION | D020 Pf8, D021 Pv4 | MalariaGEN 公开 | Pf8/Pv4 metadata 均 public_open | **不下载**，避免 >GB 级 VCF 全量 |
| STAGE/ANNOTATION/STRUCTURE | D022–D025 | GEO/PRIDE/InterPro/AFDB | 均为 public API | **不下载**，仅在锁定 orthogroup/候选后按需查询 |

---

## 5. 初始物种树可行性

- **已知框架：** 18 物种的分型与 Hamilton et al. 2017 Fig.1 谱系一致 — Laverania (Pf, Pr, Gaboni), vivax/knowlesi clade, malariae clade, ovale clade, rodent malaria, avian, piroplasms, coccidian outgroup。
- **可行性：** 12 个 FROZEN 物种均有 ≥3,600 蛋白编码基因，足以构建 **高可信单拷贝 orthogroups** 的核心物种树；P. ovale 降级种不计入核心树。
- **系统发育控制承诺（04 §3/§4）：**
  - 主分析以 **species 为统计单位**，基因不膨胀样本量；
  - 生态/宿主关联必须用 **phylogenetic comparative** 或等价控制；
  - 必做 **leave-one-clade-out**（Laverania / vivax-clade / rodent / piroplasm 各去一支）；
  - 顶质体/线粒体与核基因组分开建树与组成计算。
- **风险：** 长枝 (Toxo, piroplasms) 与 AT 极端种的 saturation；M1 将报告 gene-tree conflict 与 compositional heterogeneity 对拓扑的影响。

---

## 6. Gate G0 评审要点与证据钩子

| 证据门问题 | 最低证据 | 当前证据钩子 | 是否满足 |
|---|---|---|---|
| Are the data, species contrasts and public access sufficient to start? | ≥8 Plasmodium; ≥2 独立对照潜力；≥2 深挖谱系；无私有数据 | 12 FROZEN + 2 降级可用；≥3 对照候选；Pf + Pv/Pk/Pb 均公开可扩展；D001–D005 已校验 | **GO** |

**Gate A 预判（仅作风险提示，非结论）：**
- 2 个可靠的自然实验候选（vivax-clade vs Laverania，rodent/avian vs Laverania）很可能通过“GC 水平/区室重复性”检验；
- ovale 碎片化与 SP006 注释滞后是已知的降级触发器，若 M1 发现主要信号完全由 scaffold N50 或 annotation gene-count 差驱动，则触发 `STOP/DOWNGRADE` 至 clade-specific 设计 — 这属于**成功的项目管理**，不是失败隐瞒。

---

## 7. 独立组成对照建议（不锁定机制）

供 M1 直接检验的 **两个主对照 + 一个外群梯度**（均需在 M1 以 `GENOME / CDS / GC3 / 4D / compartment / LCR masked vs unmasked` 双轨复现）：

1. **对照 A — Laverania 高AT vs vivax/knowlesi 中GC**  
   `Pf (19.5%) + Pr (18.5%) + Gaboni (18.0%)` vs `Pv P01 (~42%) + Pk (38.5%) + PcynB (40.5%) + Pcoatneyi (39.5%)`  
   生态多样但同为灵长类寄生虫，可在控制宿主大类后检验组成差异是否为宿主生理的简单传导。

2. **对照 B — 啮齿/禽疟原虫 vs Laverania（更远缘但保留 Apicomplexa 核心）**  
   `Pberghei (22.0%) + Pyoelii (21.5%) + Pchabaudi (23.5%) + Prelictum (18.5%)` vs Laverania  
   借助不同宿主与传播生态的独立分支；*P. relictum* 的 18.5% GC 提供“非 Laverania 低GC”对照，利于区分 AT 偏好是单次起源还是多起源。

3. **外群梯度 — Piroplasm/coccidian 的 GC 高值锚点**  
   `B. bovis (41.5%), B. microti (36.5%), T. parva (34%), T. gondii (52.5%)`  
   不作为 Plasmodium 内对照的替代，仅用于检验组成-谱系祖先重建的稳定性与闭合数据方法的敏感性（Aitchison/CLR）。

> M1 必须对每个对照报告 **effect size + uncertainty + leave-one-clade-out + LCR 屏蔽前后一致性**，否则不计入 Gate A 证据。

---

## 8. 领域专家待确认清单

按 `12_COLLABORATION_AND_EXPERT_INPUT.md` 需在 G0 评审时请寄生虫学家书面确认：

1. **P. vivax 株选择：** 是否同意以 **P01 (GCA_900093555.2)** 作为 vivax 代表，Salvador-I 仅作灵敏度对照？（M0 倾向 P01：基因集更大、与 PlasmoDB-71 一致、近期文献主流。）
2. **ovale 双种的可解释性：** GH01 / wallikeri 的碎片化 scaffold 注释是否允许在 Figure 1 以“降级可用”呈现，或应直接移至附录仅作提及？
3. **第二深挖谱系优先级：** 在 *P. knowlesi* (人-猴共患 + 单细胞 + Tn-seq) 与 *P. berghei* (全生命周期单细胞 + 啮齿模型) 之间，专家更倾向哪个作为**第二**深挖支以支撑阶段/必需性交叉？
4. **P. gaboni / P. relictum 的表型可及性：** 是否有公开表型/培养/媒介数据足以支撑功能叙事，或应将其定位为“发现层仅用于组成/orthogroup 拓扑”的辅助种？
5. **M1 的可失败声明：** 若 M1 未能找到两个独立对照，是否同意按 ROADMAP 触发 `PIVOT-SCOPE` 收缩至单一类群而非追加数据抢救高目标？

---

## 9. 最小 M1 启动任务（通过 G0 后）

以下任务包均对应 `T05_TASK_PACKET.md` 单一边界交付，子 Agent 不修改 Gate：

| task_id | Aim/WP | 输入 (已冻结) | 方法边界 | 必交输出 | 验收 |
|---|---|---|---|---|---|
| M1-01 | Aim1/WP1 | 12 FROZEN genomes + VEuPathDB 71 GFF | 全基因组/区室/核-线粒体-顶质体分开的 GC, GC1/2/3, 4D, 二/单核苷酸, homopolymer, LCR 双轨 | 组成地图表 + Figure 1 审计图初稿 | N50/总长/基因数/BUSCO 核查通过；LCR 屏蔽前后方向不反转 |
| M1-02 | Aim1/WP1 | 同上 + orthogroups | 单拷贝核心 orthogroup 物种树；祖先 GC 重建；gene-tree conflict 报告 | 物种树 (newick + bootstrap) + 组成祖先状态 | leave-one-clade-out 敏感性；与文献树冲突需解释 |
| M1-03 | Aim1/WP1 | 同上 | 区室/端粒距离/重复/mappability 对组成的混杂检验 | 伪影排除报告 | 若主信号可被区室或注释完全解释则标记 STOP 信号 |
| M1-04 | Aim2预备/WP2 | 同上 | 同义密码子使用基线描述 (RSCU/ENC/PR2/PCA 仅作 QC) | QC 表 | 不作适应性断言 |

**禁止：** 在 M1 提前跑 ML 分类器宣称跨物种泛化、或拉取 MalariaGEN 全量 reads。

---

## 10. 风险、局限与降级路径

| 风险 | 影响 | 缓解（已做/待做） |
|---|---|---|
| SP009/SP010 ovale 无注释染色体 | Gate A 独立对照数虚高 | 已降级为 CANDIDATE；M1 不纳入 Gate 证据分子，仅附录 |
| SP006 P. cynomolgi 较新 Chr 未注释 (N50 3.1 Mb) | 误以碎片化惩罚真实连续性 | 保留 GCF 注释版为主；GCA 仅待重注释后晋升 |
| P. vivax P01 vs Salvador-I 株差 | 亚端粒家族数差异易被误判为物种差异 | 双株对照 + PlasmoDB 交叉检查 |
| D016 SRA raw reads >100 GB | 磁盘与时间失控 | 绝不自动下载；仅在 Gate C 明确需要 re-calling 时申请增量 |
| 长枝 + AT 极端导致的系统发育伪影 | 组成祖先重建偏倚 | 必做模型敏感性（compositional heterogeneity-aware 模型 + 非平稳模型对照） |

**降级路径（与目标期刊解耦）：**
- `A Pf-centred + cross-species 残差预测翻译应激 → Nature Microbiology`
- `B 机制弱但资源完整 → Nature Communications / PLoS Pathogens / iScience`
- `C 中性传导为主 → 比较基因组学期刊`
- `D 可塑性转向 → 若重复/CNV 证据足够仍可维持高目标`

---

## 附录 A. 文件与校验索引

```
data/metadata/species_assembly_audit.tsv          19 行  (18 物种 audit)
data/metadata/frozen_assembly_candidates.tsv      19 行  (冻结对照表)
data/metadata/taxonomy_host_ecology.tsv           19 行  (18 物种 + host 行)
data/metadata/uniprot_hemoglobin_manifest.tsv      17 行  (16 fasta)
data/metadata/nar2017_ena_accessions.tsv          285 行 (284 ENA ERS)
data/raw/ncbi-datasets/species-audit-20260903/    18 JSON + DOWNLOAD_MANIFEST.json
data/raw/ncbi-datasets/genomes/                   18 zips, 454 MB
data/raw/ncbi-datasets/<ACC>/                     18 dirs, each DOWNLOAD_MANIFEST.json + unpacked genome/gff/cds/protein
data/raw/veupathdb/PlasmoDB-71/                   46 files, 554 MB, DOWNLOAD_MANIFEST.json (46 sha256)
data/raw/ncbi_taxonomy/                          54 xml/json
data/raw/uniprot/hemoglobin/                     16 fasta
data/raw/nar2017/                                70 MB, 3 files + manifest
data/checksums/                                  6 sha256 manifests, all `sha256sum -c` OK
docs/manifests/data_registry.tsv                 D001–D005 DOWNLOADED_VERIFIED (2026-09-03)
```

**Checksum 验证命令（可复现）：**
```bash
sha256sum -c data/checksums/ncbi-genomes-20260903.sha256
sha256sum -c data/checksums/veupathdb-PlasmoDB-71.sha256
sha256sum -c data/checksums/species-audit-20260903.sha256
sha256sum -c data/checksums/taxonomy-20260903.sha256
sha256sum -c data/checksums/uniprot-hemoglobin-20260903.sha256
sha256sum -c data/checksums/nar2017-20260903.sha256
```

## 附录 B. 版本与引用

- NCBI Datasets CLI 18.36.0, UniProt release 2026_03, VEuPathDB PlasmoDB Release 71, NAR 2017 DOI 10.1093/nar/gkw1259 (PMC5389722)
- 冻结日期：2026-09-03；下次冻结变更需 `T01_CHANGE_REQUEST.md`
- 许可证：NCBI terms, VEuPathDB terms, UniProt license, PMC open article — 详见 `data_registry.tsv`

---

**状态更新建议（待 G0 评审后由 Maintainer 执行）：**
- `docs/manifests/project_state.yaml`: `status: M0_AUDIT_COMPLETE` / `last_completed_task: M0_READINESS_AND_DATA_AUDIT.md v1.0` / `next_action: G0 review → M1-01..04`
- `docs/manifests/gate_status.tsv`: `G0_READINESS` → `REVIEW_REQUESTED` (evidence_refs: 本报告 §2–§4)
- `docs/manifests/evidence_ledger.tsv` / `claim_registry.tsv`: 仅记录数据可得性与冻结事实，不新增生物学 claim

*— End of M0 Readiness and Data Audit v1.0 —*
