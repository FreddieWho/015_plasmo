# Manuscript draft v0.1 — PLASMODIUM-C2F (M5终审后手稿打包)

> 状态：DRAFT。措辞上限 = claim_registry allowed_language（CLM08/09 CANDIDATE C2，CLM10 hypothesis）；标题禁区：no repeat-first、no dual-axes、no co-option-proven、no supply-collapse。每个实质段落末尾附 claim 追溯 tag。Figure 编号见 M5_FINAL_REVIEW §3 + storyboard §11。图版文件待图包 worker 回包后挂载（manuscript/figures/）。

## Title（候选，投稿前由用户定稿）

- T1: Genome-composition transitions propagate into Plasmodium protein-sequence space and partition regulatory architectures
- T2: Composition-driven protein-sequence turnover across Plasmodium lineages is architecturally partitioned and distinct from adaptive substitutions
- T3（弱版备用）: A composition-to-protein map of Plasmodium genome evolution with a bounded regulatory layer

## Abstract（~250 words，证据顺序 per storyboard §9）

1. 领域问题：跨物种基因组组成差异长期混淆功能解释，疟原虫 AT 极端组成与其耐药/阶段适应的关系不清 [背景]。
2. 未解缺口：组成差异中多少是机械传导、多少是功能残差，缺乏系统发育控制下的定量分解。
3. 方法原则：18 物种组成地图 + 跨层反事实 + branch-aware 位点分析 + 蛋白架构分区 + 生命周期/药物双语境映射 [M1/M2/L2/M4R]。
4. 主要发现：碱基背景解释密码子差异 84.5%、氨基酸差异 78.7%（Fig2）；组成耦合位点在 informative branches 上同向替换率 0.791（OR=2.25，LOO 7/7），落在结构化保守蛋白且与 14/14 已知耐药位点分离（Fig2–3）[CLM08]；AP2 家族 Asn 富集 OR≈35（正式集复现），但家族级 lifecycle-state 后果未确立（Fig4）[CLM09]；ART 急性/慢性分区描述成立、机制未立（Fig5）[CLM10-hypothesis]。
5. 功能闭环：单基因座实证（GCN5 poly-Asn repeat，外部阳性对照）+ 群体层耐药位点分离（Pf8）[L5半]。
6. 意义：长期组成演化限定短期应激响应空间；co-option 与 ART-MOA 仍为 hypothesis，需扰动实验裁决 [L6待]。

## Results

### Fig1 独立组成转换（M1；CLM01）
- 18 物种 genome GC 18.2%（Pgaboni）→52.3%（Tg 外群）；Laverania+avian 低、vivax+piroplasm 高两组独立对照（grade A）；7/7 LOO 保留梯度；N%/contig 非混杂（|r|<0.4）。
- 失败/边界：SP009/010 无 CDS（16 物种下游）；ovale 降级；mito/api 缺失为组装收录非生物学。

### Fig2 传导 + branch-aware 桥（M2 + D1/D3；CLM02/03/08）
- GC3 解释密码子 chi2 偏离中位 84.5%（284/288 q<0.05，效应碎片化）；CLR~genome_gc 中位 R² 78.7%（14/20 斜率显著=背景本身）；7 solid + L provisional；13 robust GC-coupled AA。
- Branch-aware：informative branches 总体 0.583；coupled 0.791 vs uncoupled 0.628（Fisher OR=2.25，p=9.5e-10）；multi-nt 0.667 > 1-nt 0.546（可达性梯度）；terminal 相关保留为描述层。
- Adaptive 分离（弱版措辞）：14/14 耐药位点跨种保守、零耦合；conservation 梯度（cons=1.0 耦合 0% vs 0.5 处 9.1%）；适应性替换 Grantham 横跨全谱。不称双适应轴。

### Fig3 架构分区 + 自然实验 + 群体边界（D2b/D6/D3/Pf8；CLM08框内）
- D2b（框内限定必写）：91 harbor 基因匹配后 LCR 零超额（MW p=0.10）、Asn-贫（0.064 vs 0.119，p=3.9e-33）、domain 富集；core-181 vs 蛋白组其余的框偏差定量（length 259 vs 472 aa；Asn 0.056 vs 0.120）。结论上限：框内 churn 不依赖 repeat；框外不宣称。
- D6：T2（AP2 4/4 p=0.0625；PUF_RNA 21/27 p=0.0030）+ T3（AP2 3/4；PUF_RNA 18/25 p=0.022）组成一致；T1 注册阴性（transition polarity 过弱）；CHROM 全程 null。规则满足靠 T2+T3。
- Pf8：34/34 耐药 markers 现群体分离（如 dhfr108 ~0.95）；CNV 率 GCH1/PM2-3 0.28、MDR1 0.20。SNP-level burden deferred 明示。

### Fig4 调控富集有界使用（A；CLM09 + CLM10-H-A hypothesis）
- AP2 正式集 OR=35.3（13.1–95.0）逐值复现 regex 34.3；宽 RNA-binding 衰减至 1.67；chromatin 正式集转弱阳 1.98；sexual 集与 AP2 完全共线（不 double-count）；严格 PUF n=2 不可判；transcription 正式集退化（单 anchor，非反驳）。
- 反事实：logit 调整 OR=1.97 存活；但分层仅 high-LCR 显著（部分成立）。
- GCN5 单基因座实证（外部对照，cite-and-complement，不碰 repeat-first）+ Pb 家族 map（描述性，9 blocked orthologs Asn 高但对照 n=2 不可检验）。
- 诚实阴性框：liver 极 detection 推翻排除 / zygote null / gametocyte 弱方向（BH 未过）/ PfAP2-P DE 方向不一致 + bound 基因 Asn-贫（bounding）。**不收 shortlist，不宣称 co-option。**

### Fig5/Extended ART 分区（C；CLM10-H-C hypothesis + CLM04 boundary）
- Acute：microarray 16/27 组 stageadj 存活（Dd2 稳定、Cam3II 混合；弱 stage 参考声明 + length 衰减声明）；scRNA 独立同向（mRNA 层）；persistence 第三向（dormancy）。
- Chronic：Mok 蛋白负（5/7，最强 −0.129）vs GSE59099 转录 null → 节约在蛋白层；DiD 微弱（0.02–0.07）。
- dTE 规则2：方向一致（agree 0.924）+ q 仅 1 基因 → 弱支持保留。M4R-X01 gap 明示：无 acute translation×K13 公共数据。

## Discussion（有界）
- 主张上限：进化总原则（框内）+ 调控富集模式 + 单基因座实证 + ART 分区描述。Association 上限，无机制动词。
- 竞争解释逐条保留：length/表达/annotation 偏倚；detection 偏倚（liver）；stage 混杂残余（弱参考）；IEA 电子注释循环（L-010 已量测：curated 非IEA 层 n=10 欠定）。
- 升级/证伪条件：L-010 已回（AP2 未消失：curated 欠定，T01 维持）；家族级 perturbation 表达链（H-A 转正）；acute translation×K13 数据（H-C 转正）；D2b 框外证伪（Fig3 降级）。

## Methods（按图 pipeline + 版本 + 种子）
- 环境：~/.conda/envs/sc（pandas 2.3.3/scipy 1.11.3）+ LD_PRELOAD（D-011）；seed 20260905；MWU 一律 asymptotic（ties 组 exact 组合爆炸教训 D-037）。
- M1：k21 sketch NJ + squared-change parsimony祖先 + 7 LOO；M2：GC3 null/CLR-OLS/BH双轨/LCR双轨；L2v2：Grantham + QC + IQ-TREE ML；D1：Fitch + branch GC + Fisher；D2b：length-decile×OG-age 8层匹配；D6：mafft --auto + 预注册 pole；A：Fisher+Woolf CI+BH + logit；C：spearman + lenadj + noSchiz + BH per contrast。
- 数据：GEO（GSE75795/120448/134268/120488/225340/59099/222586/220039-Pf/151189/226632/58402）+ Zenodo（K13 scRNA 20344254；Pf8 18681980）+ EuropePMC suppl（GCN5/Mok/AP2-P/PbApiAP2/Li2024）+ Malaria Cell Atlas + PlasmoDB-71 + NCBI frozen accessions（D001–D005）。版本/sha 见 data/checksums/*.sha256 + 各 input_manifest.tsv。
- 排除项：SRA bulk raw、FASTQ、AlphaFold bulk、作者私发数据（M4R-X01、Puf1 SI 均未碰）。

## Data availability / Code availability
- 公共数据 accession 表（manuscript/tables/ 由图包生成）+ 本地复现路径（data/raw + derived + checksums）+ 环境锁（docs/manifests/environment.sc.*）。

## Claim 追溯附录（段落→claim→evidence）
- Fig1 段→CLM01→EVID-M1-001/002/003；Fig2→CLM02/03/08→EVID-M2-001/002 + EVID-M4R-D-001/002；
- Fig3→CLM08→EVID-M4R-D-002；Fig4→CLM09→EVID-M4R-A-001/002/003（含 GO 行）；
- Fig5→CLM10-H-C + CLM04→EVID-M4R-C-001；标题/摘要→不超过 CANDIDATE/hypothesis。
