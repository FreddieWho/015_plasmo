# M4R-A claim impact (ROUND1, 2026-09-05)

## CLM09（regulatory enrichment 候选）— 加强但未升级

- 富集（Asn top-decile q90=0.184，n=5285；**product-description regex 衍生基因集，非 GO/InterPro 正式集**）：ApiAP2 OR=34.3（13.8–84.9）q=1.5e-16；
  transcription_reg OR=3.84（2.41–6.12）；sexual_gametocyte OR=3.71（1.94–7.12）；
  CCR4-NOT OR=9.06（2.91–28.2，n=12）；RNA_binding_broad OR=3.08（2.21–4.31）。
  chromatin_reg n.s.（OR=1.54 q=0.15）；**严格 PUF 集仅 2 基因不可判（D 线 PUF_RNA=56 宽集 OR=3.03 口径不同，不能互换）**；proteostasis n.s.
- 反事实：Logit 调整后 is_reg OR=2.59（1.69–3.97）p=1.1e-5（控制 loglen + lcr_frac）→
  富集不被长度/LCR 完全解释。但分层后仅在 high-LCR 层显著（median 0.206 vs 0.173，p=4.8e-7），
  zero-LCR 层 n.s.（p=0.16）→ “given comparable IDR” 反事实部分成立，非普遍。
- Asn-rich ≠ poly-Asn：top10 中 57% 含 polyN≥10，43% 为分散型 → 后续必须分解 tract 类型。

## CLM10-H-A（lifecycle co-option hypothesis）— 初步方向性证据，仍为 hypothesis

- GSE75795（n=1/性别，描述性）：spearman(Asn, F-vs-M logFC)=−0.088 p=5.3e-10；
  Asn-top10 中位 −0.134 vs rest +0.124（MW p=0.037，BH q=0.13 n.s.）→ Asn 富集基因向 male 偏；
  sexual_gametocyte 集向 female 偏（+0.70，p=0.035 q=0.13）。效应小、未过 BH、多重比较下保留探索标签。
- AP2-G ChIP（单研究）：bound 基因轻微 Asn 富集（MW p=0.045），bound × Asn-top10 OR=1.55 p=2.6e-4，
  bound × ApiAP2 OR=3.15 p=0.0058；AP2-I bound 反而 Asn 贫（p=6.9e-4）→ occupancy 锚弱阳性，
  非独立复现。
- GCN5 2026（cite-and-complement）：Fig2 source data 证实条件性 GCN5 扰动致生长崩溃
  （day6 −RAP ~34–37 vs +RAP ~2.0–2.3，~15×）；Fig6 为 H3 相关肽段 HAT 读数；repeat-deletion +
  Py-complement 细节在正文/Figs3–9。本轮 = 阳性对照在手，不碰 repeat-first 措辞（NT-1）。

## 禁止外推声明

- enrichment ≠ mechanism；未宣布 lifecycle co-option 成立；未提湿实验设计。
- 候选 shortlist 本轮不收（统计成立标准未达：lifecycle 端 BH 未过 + n=1 + ChIP 单研究）。

## Verdict：CONDITIONAL（A 作为 functional pillar 条件性成立，详见 M4R-A_ROUND1_REPORT.md）

## ROUND2 追加（2026-09-05）
- CLM09：AP2 正式集（IPR001471，运行时从数据解析，bg=24 完全特异）OR=35.3 vs regex 34.3，逐值复现；严格 PUF n=2 不可判（PUF1 富集/PUF2 不富集），repeatbroad null；RNA_binding 衰减但存活（1.67）；chromatin 正式集转阳（1.98）；CCR4-NOT 衰减至 n.s.；transcription 正式集退化；sexual_gametocyte formal == AP2 集（不 double-count）；logit 调整 OR=1.97。anchor 发现用 product 文本（Asn 检验独立于选择），正式性弱于人工 GO——已声明。
- CLM10-H-A：三极中 liver BH-pass 但 detection 对照推翻（NoGFP 检出不对称 p=3.7e-11；检出控制后 n=80 null）→ 排除；zygote null；gametocyte 弱方向维持探索。state-consequence 标准未达，不收 shortlist。
- Perturbation：PfAP2-P 身份已验；DE 方向不一致；bound 基因 Asn-贫（p~1e-35/1e-16，top10 OR≈1）→ bounding 证据（家族组成 vs 单因子 regulon 分层）；PbApiAP2 map 落盘（27 尝试/11 verified/9 transmission-blocked），交叉 n=2 不可检验。
- Verdict：CONDITIONAL 不变。详见 M4R-A_ROUND2_REPORT.md。
