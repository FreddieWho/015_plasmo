# Sub-Agent Task Packet — MANUSCRIPT FIGURES D-side (Fig1–3)

- task_id: MS-FIG-D
- parent milestone: M5终审后手稿打包（D-039 Figure 映射 §3；storyboard §10/§11）
- linked claims: CLM01/02/03（Fig1–2）+ CLM08（Fig3，CANDIDATE C2 上限，框内限定）
- job: 只绘图 + 图例草稿，不做新分析。从现有结果 TSV 画 Fig1–3 每个 panel（PNG 300dpi），每图配生成脚本 + panel 数据表 + legend 草稿 TXT。
- inputs (ALL IN HAND):
  - Fig1: data/derived/WP1/M1-01_composition/M1-01_composition_table.tsv + M1-01_compartment_breakdown.tsv；M1-02_species_tree.nwk + M1-02_ancestral_gc.tsv + M1-02_independence_grade.tsv + M1-02_leave_one_clade_out.tsv
  - Fig2: data/derived/WP2/M2-01_codon/M2-01_background_explained.tsv（84.5%）+ M2-02_aa/M2-02_background_explained.tsv（78.7%）；M4RD_D1_concordance_summary.tsv（coupled 0.791 vs uncoupled，CI 列齐）+ M4RD_D1_LOO.tsv（OR 1.71–4.35 范围带）；M4RD_D3_multidim.tsv + M4RD_D3_matched_controls.tsv（14/14 uncoupled + conservation 梯度）
  - Fig3: M4RD_D2b_matched.tsv（harbor vs 匹配背景 LCR/Asn）+ M4RD_D2b_framebias.tsv（采样框偏差定量）+ M4RD_D2b_verdict.txt（SURVIVES）；M4RD_D6formal_summary.tsv + M4RD_D6formal_pertransition.tsv（T2/T3 + T1 注册阴性 + CHROM null）；M4RD_D4_Pf8boundary.tsv（markers 分离）+ M4RD_D4_Pf8cnv.tsv（CNV 率）
- panels:
  - Fig1: (a) 18 物种 genome GC bar + 区室分层；(b) 物种树 + 祖先 GC + independence grade 标注；(c) LOO 稳健小注
  - Fig2: (a) GC3→codon 84.5% / genome_gc→AA 78.7% 传导条；(b) branch-aware coupled-vs-uncoupled forest（OR=2.25 + LOO 范围）；(c) 耐药位点分离示意（14/14 uncoupled + conservation 梯度，弱措辞）
  - Fig3: (a) D2b harbor-vs-matched（LCR 零超额 + Asn-贫 + 框偏差注）；(b) D6 三 transitions（T2/T3 一致 + T1 阴性 + CHROM null）；(c) Pf8 群体层（markers 分离率 + CNV 率）
- legend 措辞上限：CANDIDATE/C2 内；禁 repeat-first/dual-axes/co-option-proven/supply-collapse；Fig3 必带框内限定句；T1 阴性必须出现
- QA §7 per figure: 一句话主张 + 数据表 + 生成脚本 + N/统计单位 + 效应/区间 + 敏感性/Extended 指向 + 非颜色依赖标注 + 物种/基因ID/版本
- outputs: manuscript/figures/Fig{1,2,3}{a,b,c}.png + manuscript/figures/scripts/mk_fig{1,2,3}.py + manuscript/figures/panels/*.tsv + manuscript/figures/legends/Fig{1,2,3}_legend.txt + FIG_D_REPORT.md（含每 panel 输入 sha + 与 claim 对应关系）
- guardrails: 超时套壳；只读输入不写回 derived；matplotlib 无中文字体用英文标注；色板 colorblind-safe + 形状双编码
- return: 输出路径 + 每图一句话主张 + 与 claim 对应 + 失败/限制
