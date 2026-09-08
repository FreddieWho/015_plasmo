# Sub-Agent Task Packet — MANUSCRIPT FIGURES A/C-side (Fig4–5 + Extended + Tables)

- task_id: MS-FIG-AC
- parent milestone: M5终审后手稿打包（D-039 Figure 映射 §3；storyboard §11：Fig4 有界 A，Fig5/Extended C）
- linked claims: CLM09（CANDIDATE；AP2 强/宽 RNA-binding 弱/chromatin 弱阳/PUF 未定）+ CLM10-H-A/H-C（hypothesis）+ CLM04（boundary context）
- job: 只绘图 + 图例草稿 + 主表索引，不做新分析。
- inputs (ALL IN HAND):
  - Fig4: M4RA_interpro_enrichment.tsv（正式集 vs regex 并排：AP2 35.3/34.3、chromatin 翻转、sexual==AP2 共线、transcription 退化）；M4RA_counterfactual_logit.tsv + stratified（logit 1.97 + high-LCR-only）；M4RA_gcn5_fig2_snapshot.tsv（~15× 生长缺陷）+ inventory；M4RA_pbapiap2_family_map.tsv（9 blocked orthologs 描述性）；M4RA_lifecycle_tripole.tsv + M4RA_liver_detection_aux.tsv（三极阴性报告框）；**M4RA_l010_enrichment.tsv（如运行时存在则加 GO 复核行，否则 legend 注 pending L-010）**
  - Fig5: M4RC_acute_chronic_partition.tsv（asn_frac 行：acute 正两端 / Mok 负 / GSE59099 null / DiD 微弱）；M4RC_dTE_upgrade.json（agree 0.924 + q 仅 1 基因，规则2）；M4R-X01 gap 声明文本（M4R-C round-1 report 有）
- panels:
  - Fig4: (a) 正式-vs-regex 富集 forest（OR+CI；sexual 行标注共线不 double-count；PUF 行标注 n=2 不可判）；(b) 反事实分层（high-LCR-only 成立声明）；(c) GCN5 单基因座实证 + Pb map 描述 + 三极阴性诚实框（liver 排除/zygote null/gametocyte 探索）
  - Fig5: (a) acute-vs-chronic 分区 forest（16/27 存活注 + 弱 stage 参考声明 + length 衰减声明）；(b) dTE 规则2盒 + M4R-X01 gap 盒
- legend 措辞上限：enrichment≠mechanism；不宣称 co-option；mobilization/conservation 只作 hypothesis；C 为第二语境不独立承重
- QA §7 per figure（同 D 包）
- outputs: manuscript/figures/Fig{4,5}{a,b,c}.png + scripts/mk_fig{4,5}.py + panels/*.tsv + legends/Fig{4,5}_legend.txt + manuscript/tables/TableS{1..6}.tsv（D1/D2b/D6/interpro/tripole/C-partition 现表复制 + tables/README_index.md 注每表来源 sha）+ FIG_AC_REPORT.md
- guardrails: 同 D 包；liver 面板必须同时画 detection 对照结论（排除章）；GSE 及 GEO accession 号标注版本
- return: 输出路径 + 每图一句话主张 + 与 claim 对应 + 失败/限制
