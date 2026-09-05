# M4R-A ROUND1 — README

- task: M4R-A-ROUND1（Route A functional/lifecycle pillar，第一轮信息增益）
- date: 2026-09-05； env: `~/.conda/envs/sc` python + LD_PRELOAD libstdc++（D-011）
- script: run_m4ra.py（params.yaml；input_manifest.tsv；checksums.sha256）

## 输入与版本

- PlasmoDB-71 Pf3D7 GFF + AnnotatedProteins.fasta（最长 isoform/基因，n=5285）
- L1d 基线：L1d_keyword_enrichment.tsv（关键词版起点，非终点）
- GSE75795 RAW.tar（134MB）：PlasmoDB 9.1 raw bedGraphs（plus+minus）× 3 样本（male GFP / female mCherry / male double-positive）；9.1 坐标，v71 基因区间求和；每性别 n=1 → 描述性
- ChIP BED（rep-consensus overlap）：GSE120448 AP2-G S/R/G；GSE134268 NCC/SCC；GSE120488 AP2-G/AP2-I schizont；基因映射 = body 或 TSS 上游 ≤1.5kb；v71 坐标近似
- GCN5 2026 suppl：MOESM3（Mascot 蛋白组 715 行）/ MOESM4（primers 31 行）/ MOESM7（Fig2–9 + S6 source data）
- 缺：GSE222586（mosquito midgut）/ GSE220039（liver）本地未入库 → 未算，已报缺；PfAP2-P/PfPuf1KO/Hsp110c/PbApiAP2 processed 表未入库 → 未算

## 方法要点

- 蛋白特征：Asn fraction、poly-Asn 最大 tract、LCR proxy（window-20 Shannon entropy<1.8 占比，heuristic 非 InterPro）、长度
- 基因集：product-description regex（GFF 无 GO 条目，明确为 description-derived）；Fisher OR + Woolf 95%CI + BH q
- 反事实：Logit（is_reg + loglen + lcr_frac）+ LCR 三层分层 MW
- 表达：bedGraph 基因求和 → CPM → log2FC(F-vs-M)；spearman + MW（探索性 p）
- 输出表：M4RA_formal_enrichment / counterfactual_logit / counterfactual_stratified / gametocyte_counts / lifecycle_mapping / chip_consensus / chip_overlap / gcn5_inventory / gcn5_fig2_snapshot

## 局限

1. LCR proxy 是启发式 entropy 窗口，不是 InterPro/SEG；IDR 无直接注释 → 反事实是 LCR-adjusted 不是 IDR-matched。
2. Gametocyte 每条件 n=1，无 replicate → 方向性描述，不做 discovery 宣称。
3. ChIP trio = 单研究（PMID 32198457），intra-study orthogonal only；9.1→v71 坐标近似。
4. PUF 严格集仅 2 基因（描述标注稀疏），family-wide PUF 结论不可下；宽 RNA-binding 集作补充。
5. GCN5 全机制解析需论文正文；本轮只做 inventory + 生长表型快照。
