# M4R-A ROUND2 REPORT — InterPro formal sets + strict PUF + lifecycle tripole + perturbation crosswalk (2026-09-05)

- task_id: M4R-A-ROUND2； packet: docs/tasks/M4R-A_ROUND2_PACKET.md
- verdict: **CONDITIONAL（等级不变；富集臂加强，state-consequence 臂削弱；不收 shortlist；H-A 仍为 hypothesis）**
- outputs: /home/huyudi/015_plasmo/data/derived/WP4R/M4R-A/（新增 M4RA_interpro_map.tsv, M4RA_interpro_sets.tsv, M4RA_interpro_enrichment.tsv, M4RA_strictPUF.tsv, M4RA_interpro_logit.tsv, M4RA_lifecycle_tripole.tsv, M4RA_liver_detection_aux.tsv, M4RA_perturbation_xwalk.tsv, M4RA_pbapiap2_inventory.tsv, M4RA_pbapiap2_family_map.tsv, M4RA_pbapiap2_xwalk.tsv, run_m4ra_r2s1/s2/s3/s4/s5.py, 本报告；claim_impact.md 追加段；input_manifest.tsv 追加行；checksums_round2.sha256）
- method: Fisher OR+Woolf CI+BH（与 round-1 同口径）；Logit length/LCR 调整；MW 一律 asymptotic；Asn top-decile = q90 = 0.1837（universe n=5285，与 round-1 一致）

## 1. 输入版本
- NCBI GCF_000002765.6 genomic.gff（sha 76206d1e…，12.5MB）：已验 == Pf3D7（gene 行 locus_tag/Name=PF3D7_*；CDS 行 Dbxref 含 InterPro:IPR* + locus_tag）。映射：CDS locus_tag→基因，InterPro 取并集。5282 基因有 CDS，其中 3622 含 ≥1 IPR。
- GSE222586 zygote 矩阵+metadata（sha 74a98bb9…/7231024a…）；GSE220039 Pf 肝期 counts（sha 3abc4f11…）。
- PMC10627835 MOESM3（sha cd804370…）/MOESM5（sha b117c636…）：sheet 明示 pfap2-p truncated / PfAP2-P bound → **身份已验，可用**。
- PbApiAP2 PMC5241200_supp.zip（sha 5197394f…，28MB）：mmc1/7 pdf 未解析（无 pdf 工具，map 备注）；mmc2–4/6 xlsx 已 inventory。
- Puf1 SI：按 packet 未重试（D-038 不可达定论维持）。

## 2. IPR 号全部运行时从数据解析（未背号码）
- AP2 domain = **IPR001471**（24/24 product-AP2 基因携带；全基因组 bg=24 → 完全特异）。正式 AP2 集 n=24。
- Pumilio = **IPR001313**（2/3 product-PUM 基因携带；bg=2）。严格 PUF 集 n=2。
- 其余正式集 = product-text 类内 per-IPR Fisher（BH q<0.05 且 OR>2 且 n≥3）选出的 anchor IPRs（方法见 run_m4ra_r2s2.py）。**方法声明：anchor 发现用 product 文本，Asn 检验独立于选择步骤；因此这是 InterPro-anchored 定义，不是独立 GO 注释——正式性弱于 GO/InterPro 人工注释，强于纯 regex。**

## 3. 正式富集 vs regex 基线（并排）
| 集 | 正式 OR (95%CI) q | regex 基线 OR q | 判定 |
|---|---|---|---|
| ApiAP2 (IPR001471, n=24) | 35.3 (13.1–95.0) 1.2e-14 | 34.3 (13.8–84.9) 1.5e-16 | **复现，几乎逐值一致** |
| 严格 PUF (IPR001313, n=2) | 9.0 (0.56–144) n.s. | 9.0 n.s. | **不可判（集太小）；CLM09 不得写 PUF** |
| PUF_repeatbroad (79) | 1.16 (0.57–2.33) n.s. | — | 宽 repeat 集不复制；通用 repeat IPR 非 PUF 证据 |
| RNA_binding_broad (467, 23 anchors) | 1.67 (1.27–2.19) 9.4e-04 | 3.08 (2.21–4.31) 8.3e-09 | 方向一致、显著但衰减 |
| chromatin_reg (176, 17 anchors) | 1.98 (1.33–2.94) 3.8e-03 | 1.54 (0.94–2.53) n.s. | **翻转：正式集显著，regex 不显著** |
| CCR4_NOT (8, 2 anchors) | 3.0 (0.60–14.9) n.s. | 9.06 (2.9–28.2) 1.2e-03 | 衰减至 n.s.（n 小） |
| transcription_reg (4, 1 anchor) | n=4, 0 in top | 3.84 significant | **正式集退化（单 anchor），不是反驳** |
| sexual_gametocyte formal | == AP2 集（唯一 anchor 即 IPR001471） | 3.71 | **与 AP2 完全共线，不许 double-count** |
| proteostasis (318, 34 anchors) | 1.12 (0.78–1.61) n.s. | 1.17 n.s. | null 维持 |
- Logit（is_reg_ip=AP2∪PUFstrict∪chromatin∪CCR4NOT，控制 loglen+lcr_frac）：调整 OR=**1.97 (1.32–2.95)** p=9.4e-04（round-1 regex 版 2.59）。注：lcr_frac 系数极大（分离迹象），is_reg_ip 估计不受影响但报告备查。

## 4. 严格 PUF 表（M4RA_strictPUF.tsv）
- PF3D7_0518700 PUF1：Asn 0.272，polyN 33，top10=True。
- PF3D7_0417100 PUF2：Asn 0.109，polyN 2，top10=False。
- n=2 → 不可判。CLM09 禁写 PUF（维持 D-037 措辞）。

## 5. Lifecycle 三极（BH 统一；M4RA_lifecycle_tripole.tsv）
- gametocyte（n=4949 表达基因）：spearman −0.088 q=2.7e-09（Asn 富集向 male 偏）；MW top10 q=0.081 n.s.；AP2 formal MW q=0.13 n.s.（中位 −0.954 vs +0.096，方向与 Asn 一致）。
- zygote（T0/T2/T4/8/12/20h，125 细胞伪bulk）：全部 n.s.（最优 q=0.073）。null。
- liver（GFP-vs-NoGFP，Day2/4/6 + Day0）：spearman +0.10 q=4e-13、MW q=3e-10、Day4/Day6/Day0naive 全部 BH 通过；AP2 formal Day0naive q=0.034，Day4 q=0.074。**但 detection 对照推翻**：NoGFP 中 Asn-top10 检出率 7.3% vs rest 3.4%（MW p=3.7e-11）；双条件各 ≥3 样本检出基因仅 84/5720，其上信号消失（spearman −0.071 p=0.53）。→ **liver 极判定为 detection 驱动，排除出 state-consequence 证据。**
- state-consequence 标准（≥2 极同向 + ≥1 BH 通过）：gametocyte（弱方向）+ zygote（null）+ liver（排除）→ **未达。不收 shortlist。**

## 6. Perturbation 交叉（M4RA_perturbation_xwalk.tsv）
- PfAP2-P truncation DE-vs-Asn：16h +0.10 / 40h −0.05 / 8h −0.11 / 30h −0.03 → **方向不一致，无 coherent 扰动→Asn 关系**。
- **PfAP2-P bound 基因 Asn-贫**：16h 中位 0.092 vs 0.122（p=7.1e-35）；40h 0.090 vs 0.120（p=9.1e-16）；bound×Asn-top10 Fisher OR≈1.0 n.s. → AP2-P 调控子与 Asn 轴分离（bound 甲方），与 AP2-G ChIP 弱阳性形成因子间对照。bound×AP2 formal OR=4.6（16h）→ 自调控，合理。
- PbApiAP2 家族 map（M4RA_pbapiap2_family_map.tsv）：27 KO 尝试，11 KO-verified 具 Pf ortholog，其中 9 transmission-blocked；blocked ortholog Asn 中位 0.229（基因组 ~0.118）但对照仅 n=2 → **不可检验，仅描述记录，不作证据**。Data_S2/S3/S5（DE/共表达/双 KO 组学）已 inventory 未解析（Pb 表达≠Pf 证据，按 packet 止步于 map）。

## 7. Verdict 理由（CONDITIONAL 不变）
- 加强：AP2 正式集逐值复现（35.3 vs 34.3）；chromatin 正式集转阳；logit 调整存活（1.97）；PfAP2-P 身份验后 bounding 证据（bound 基因 Asn-贫）+ Pb 家族 map 落盘。
- 削弱：严格 PUF 仍 n=2；repeatbroad null；liver 极 detection 推翻；PfAP2-P DE 方向不一致；state-consequence 标准未达。
- 净：CLM09 维持候选（AP2 强 + 宽 RNA-binding 弱 + chromatin 新弱阳 + PUF 未定）；CLM10-H-A 维持 hypothesis；shortlist 空。

## 8. 什么能改变 verdict
- →STRONG：第二极 BH-pass state consequence（非 detection 驱动）+ 同背景 perturbation 表达链 + 严格 PUF 集可判。
- →WEAK/BOUNDARY：真 GO 人工注释下 AP2 富集消失（本轮 anchor 法不能排除此可能，已声明）。

## 9. QC / 失败记录
- sc 环境 + LD_PRELOAD（已知坑）；`M["prod"]` 列名冲突（round-1 `D.product` 同款 bug 第二次出现——建议入全局教训：DataFrame 描述列禁用 product/count/median 等保留名）；padj 列字符串混入→to_numeric；mmc pdf 未解析（无工具，如实备注）；GSE220039 12 万行分块读（20000/chunk）一次通过；GSE222586 `;`分隔；单步均一次通过，无 UNRESOLVED（除 packet 预设的 Puf1/AlphaFold/SRA/Pv 未碰项）。
- 未下载任何数据；身份未验数据零引用（MOESM3/5 已验才用）；Pb 数据未作 Pf 表达证据。

## 10. 禁止外推声明
- enrichment ≠ mechanism；未宣布 co-option；sexual_gametocyte formal 行不作独立证据；bound-Asn-贫是 bounding 不是反驳 AP2 富集（家族组成 vs 单个因子 regulon 是两个层次）；未提湿实验。
