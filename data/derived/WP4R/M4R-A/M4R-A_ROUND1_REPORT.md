# M4R-A ROUND1 REPORT — regulatory LC-IDR × lifecycle first-gain analysis (2026-09-05)

- task_id: M4R-A-ROUND1； packet: docs/tasks/M4R-A_TASK_PACKET.md
- verdict: **CONDITIONAL**（A 作为 functional/lifecycle pillar 条件性成立；未达 STRONG）
- outputs: /home/huyudi/015_plasmo/data/derived/WP4R/M4R-A/（README.md, params.yaml,
  input_manifest.tsv, checksums.sha256, run_m4ra.py, M4RA_*.tsv ×9, claim_impact.md）

## 1. 输入版本

- PlasmoDB-71 Pf3D7 GFF + AnnotatedProteins.fasta（最长 isoform，n=5285 基因）
- L1d 基线（keyword 版起点）；GSE75795 RAW.tar sha cdfcc1e0…（9.1 bedGraphs，求和到 v71 基因区间）
- ChIP BED：GSE120448（AP2-G S/R/G）/ GSE134268（NCC/SCC）/ GSE120488（AP2-G + AP2-I schizont），rep-consensus
- GCN5 2026 suppl：MOESM3（Mascot 715 行）/ MOESM4（primers）/ MOESM7（Fig2–9 + S6 source data）
- 缺（报缺未猜）：GSE222586 / GSE220039 未入库；PfAP2-P / PfPuf1KO / Hsp110c / PbApiAP2 processed 未入库

## 2. 关键效应量

### (1) 富集（Asn top-decile，q90=0.184；Fisher OR + 95%CI + BH q；**基因集为 product-description regex 衍生，非 GO/InterPro 正式集**——round-2 硬项）

| set | n | in_top10 | OR | 95%CI | q |
|---|---|---|---|---|---|
| ApiAP2 | 28 | 22 | 34.3 | 13.8–84.9 | 1.5e-16 |
| transcription_reg | 89 | 26 | 3.84 | 2.41–6.12 | 9.4e-07 |
| sexual_gametocyte | 45 | 13 | 3.71 | 1.94–7.12 | 8.0e-04 |
| CCR4-NOT | 12 | 6 | 9.06 | 2.91–28.2 | 1.1e-03 |
| RNA_binding_broad | 201 | 49 | 3.08 | 2.21–4.31 | 9.2e-09 |
| chromatin_reg | 131 | 19 | 1.54 | 0.94–2.53 | 0.15 n.s. |
| PUF（严格集，仅 2 基因） | 2 | 1 | — | — | n.s.（集太小不可判）；D 线 PUF_RNA=56 基因为宽集（含一般 RNA-binding），两者口径不同 |
| proteostasis | 166 | 19 | 1.17 | 0.72–1.90 | n.s. |

（基因集为 product-description regex 衍生；GFF 无 GO 条目。）

### (2) IDR/LCR 反事实

- Logit（is_reg + loglen + lcr_frac）：is_reg 调整 OR=2.59（1.69–3.97），p=1.1e-05 →
  富集不被长度/LCR 完全解释。
- 分层 MW（regulator vs rest，Asn frac）：zero-LCR 层 0.108 vs 0.099 p=0.16 n.s.；
  lowpos 层 p=0.68 n.s.；high-LCR 层 0.206 vs 0.173 p=4.8e-07 →
  超额 Asn 集中在已有 LCR 的蛋白，“given comparable IDR” 反事实部分成立、非普遍。
- Asn-rich ≠ poly-Asn：top10 中 57% 含 polyN≥10，43% 为分散型 Asn 富集。

### (3) 生命周期映射（GSE75795，4949 表达基因，每性别 n=1 → 描述性）

- spearman(Asn frac, female-vs-male logFC) = −0.088，p=5.3e-10（Asn 富集向 male 偏）。
- Asn-top10 中位 logFC −0.134 vs rest +0.124（MW p=0.037；BH q=0.13 未过）。
- sexual_gametocyte 集向 female 偏（中位 +0.70，p=0.035，q=0.13 未过）；ApiAP2 向 male 偏（−0.81，p=0.074）。
- 结论：方向性信号存在但小、BH 未过 → 探索标签，不做 discovery 宣称。

### (4) AP2-G ChIP（单研究 PMID 32198457，intra-study orthogonal only）

- 共识峰：S 195 / R 161 / G 404 / NCC 457 / SCC 277 / schizont AP2-G 214 / AP2-I 184；
  bound 基因：union AP2-G 系 799，AP2-I 208。
- AP2-G-bound 基因轻微 Asn 富集（中位 0.1189 vs 0.1181，MW p=0.045）；
  bound × Asn-top10 OR=1.55 p=2.6e-04；bound × ApiAP2 OR=3.15 p=0.0058。
- AP2-I-bound 反而 Asn 贫（0.101 vs 0.119，p=6.9e-04）→ occupancy 锚弱阳性。

### (5) GCN5 2026 inventory（cite-and-complement，不碰 repeat-first）

- MOESM3：Mascot 蛋白组（protein accession + spectral counts，715 行）。
- MOESM4：primers（31 行）。
- MOESM7：Fig2 生长曲线 source（条件性 GCN5 扰动：day6 −RAP ~34–37 vs +RAP ~2.0–2.3，
  ~15× 生长缺陷）；Fig6：组蛋白肽段 HAT 读数（126×93，H3_3_8 等修饰比例）；
  Fig9（117×18）及 Figs3–5/7–8 为 repeat-deletion + Py-complement 各面板 source data。
- 证明力陈述：单基因座 perturbation + natural-sequence（Py 短 repeat 互补）实验证据确立
  poly-Asn repeat 的功能后果（lifecycle-state consequence 按论文正文）；
  本项目引用为阳性对照 + NT-1 威胁规避（family-wide AP2/PUF axis 仍是本项目待证候选）。

## 3. Verdict 理由（CONDITIONAL）

- 强侧：AP2 OR=34（CI 下限 13.8）+ 多调控集同向 + 调整 OR=2.59 存活 + GCN5 阳性对照在手。
- 弱侧：lifecycle 端 BH 未过且 n=1；ChIP 单研究弱阳性；严格 PUF 集不可判；
  反事实仅 high-LCR 层成立；候选 shortlist 标准未达（本轮不收 1–3 候选）。
- 因此：A 是 functional pillar 的条件性候选 —— 富集真实且不完全由长度/LCR 解释，
  但 lifecycle-state consequence 尚未超出“regulators are Asn-rich”。

## 4. 什么能改变 verdict

- → STRONG：GSE222586/GSE220039 入库后阶段映射 BH 通过 + 同背景 matched perturbation
  表达数据（GCN5 位点除外）+ 正式 GO/InterPro 富集复现 + AP2/PUF 独立 occupancy 数据。
- → WEAK/BOUNDARY：InterPro 正式集下富集消失，或 lifecycle 三集方向不一致，
  或 high-LCR 层效应被 annotation/detection 偏倚解释。

## 5. QC / 失败记录

- scipy 系统 libstdc++ 缺 GLIBCXX_3.4.29 → 切 sc 环境 + LD_PRELOAD（D-011 同案）。
- pandas `D.product` 方法名冲突 → 改 D["product"]；qcut 在零膨胀 LCR 上崩 → 改 zero/lowpos/high 三层。
- bedGraph 成员名匹配 ".9.1." 误写（实为 "_9.1."）→ 修正后 6 文件基因求和成功。
- 未下载任何数据；未改 docs/manifests / DECISIONS。
