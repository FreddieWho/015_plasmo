# 项目指导文档（v1）：最初目标、路线转向、现状漏洞、升级路线与最小湿实验验证

> 写给谁：项目负责人 + 有生物学训练的专家，用来决定下一步做什么、不做什么、如果做湿实验怎么做最省。
> 不是成文草稿：这是 map，不是 manuscript。所有实质结论沿用 claim 注册表封顶（CLM08/09 候选、CLM10 假设），标题禁区不变。
> 图版：与叙事版 HTML 共用 manuscript/figures/ 的 15 张图，附在文末。

## 0. 阅读指南（一页）

- **本文档回答五个问题**：我们最初想研究什么（§1）→ 什么结果让我们换了重点（§2）→ 现在手里到底是什么证据、什么等级（§3）→ 还缺什么、能怎么升级（§4 计算 / §5 湿实验）→ 我们的创新点在哪（§6）。
- **两个比喻贯穿全文**（来自项目内部的诚实记法）：静态密码子偏好是**脆弱性地图**（谁容易被影响），动态 tRNA 供应是**天气**（实际什么时候崩）。地图可以画得很细，但不代表那天下雨——桥接两者需要实证。
- **证据等级纪律**：association（相关/共现，不谈因果）→ mechanism（扰动实验证明的因果）。当前所有新结论停在 association；只有引用别人的扰动实验才谈功能。

## 1. 最初希望研究什么

### 1.1 用户最关心的问题（PLAN.md §1–§2 原样立场，未改写）

疟原虫不同分支之间，基因组**碱基组成（GC/AT 分布）发生的转换**，哪些只是突变、历史和蛋白约束的被动后果，哪些被编码层**放大成了翻译应激程序或基因组可塑性**，并进一步影响耐药、增殖和生命周期适应？

最初的工作假设（Gate C/D 前不得当结论）：**不同谱系的组成转换被选择利用，编码了不同的翻译应激响应程序。**

五条假设的原始判据：
- **H1 多状态组成演化**：不是单一 AT 化轨迹，而是多个分支、多个区室各自的状态。
- **H2 跨层机械传导但非完全解释**：碱基背景解释大部分密码子差异和部分氨基酸差异，只剩少数稳定残差。
- **H3 序列编码的翻译应激程序**：稳定残差集中在特定生命阶段与应激模块，并与 tRNA/翻译供需、RNA—蛋白输出匹配。
- **H4 条件性可转化脆弱性**：组成演化造成的翻译依赖形成不可替代、较难逃逸的弱点。
- **H5 基因组可塑性备选轴**：翻译残差弱时，AT-rich 重复可能通过 indel/CNV/抗原变异产生价值。

### 1.2 为什么这个问题值得做（PLAN.md §6，未改写）

只有连到至少一个接口才有 Nature Microbiology 级意义：青蒿素耐受与蛋白稳态、跨物种功能可迁移性、性发育/孢子体/肝期翻译切换、可塑性导致的逃逸、寄生虫特异翻译机制。否则主动降级期刊，不靠加数据抢救叙事。

## 2. 什么结果让我们更换了重点（按时间线的诚实链）

### 2.1 H1/H2 通过，但判据被收紧

- **组成地图（Gate A GO）**：18 物种 GC 18.2%→52.3%，Laverania+avian 低、vivax 支+梨形虫支高，两组独立对照 grade A，7/7 LOO 保留梯度，组装指标非混杂（|r|<0.4）。
- **传导定量（Gate B GO，窄化）**：GC3 解释密码子偏离中位 84.5%；genome GC 解释氨基酸 CLR 变异中位 R² 78.7%。**但全基因组适应性主张在此降级**——背景斜率显著的是"背景本身"，没有第二层证据，不谈适应。
- 到这一步，H1/H2 通过，但项目主动把自己锁进 C2 天花板。

### 2.2 我们亲手杀死了 H3 跨物种版（Gate C PIVOT，负结果永久保留）

- **M3-03 干净阴性**：DHA×K13 时间序列 0 早期命中、效应 1e-4、TE 一致性 0/8；K13 mRNA 上调证明 assay 灵敏（不是死数据）。
- **I5 双轨补跑**：LCR 遮蔽 vs 未遮蔽一致性 0.9991、双轨零命中——PIVOT 对序列处理稳健，不是统计功效问题。
- 按预注册规则自动 PIVOT：跨物种纲领性主张放弃；CLM04 收窄为"Pf 特异候选机制"（SUPPORTED_CANDIDATE，跨物种语言禁用）。
- **此时已起草 T01 期刊降级**（Nat Commun/PLoS Pathog/iScience 备选）——文件保留，未删除。

### 2.3 新发现触发 M4R 有限重启（用户批准，不是翻盘）

Gate C 之后出现了三件不在旧判据里的事，构成有限重启的正当性（D-036）：
1. **L2 v2 位点层结果**：组成耦合位点反而更激进（Grantham 77.3 vs 66.9，p=0.0005），13/13 已知耐药位点跨种保守、零耦合——"组成轴 ≠ 适应轴"是新的独立角度。
2. **Asn 调控层拆解**：AAT 信号本质是**蛋白层 Asn 含量**（AP2 Asn 24% vs 基因组 12%，p=2e-12；AAT 同义选择比例无差异，p=0.73）——氨基酸组成框架接手，密码子框架降格。
3. **急性/慢性区制对**：急性饥饿动员（rho=+0.167）vs 慢性耐药节约（rho=−0.129）——同一根轴的两个区制。

新问题（重写后的中心问题，不是旧问题的翻版）：基因组组成的长期漂移，如何传导到氨基酸和蛋白序列空间？落在蛋白的什么部位？与经典适应性位点有何不同？其中是否有被生命周期或应激系统拿来用的？——**最后一个问号是开放的，答案允许是"没有"**。

## 3. 当前研究手里到底是什么（证据分级，不虚报）

### 3.1 已经成立（association 上限内，能写进论文）

- **组成-蛋白传导的进化总原则（框内）**：组成耦合位点在信息充分的进化分支上 79% 替换与该分支 GC 漂移同向（对照 63%，OR=2.25），LOO 7/7 全显著；周转落在结构化保守蛋白里（框内，采样框偏差已量化）。
- **适应轴分离（弱版）**：14/14 已知耐药氨基酸位置跨种保守、零耦合；Pf8 当代人群 34/34 markers 仍在分离——两条轴在进化上分离。
- **AP2 Asn 富集模式**：OR=35 双方法逐值复现 + GO 全集同向（IEA 驱动）。

### 3.2 只到候选/假设（能写，但必须带帽子）

- **CLM08（组成驱动分区）候选**：branch-aware + LOO + D2b + D6 + Pf8 已齐，但**框内限定**是硬条件。
- **CLM09（调控富集）候选**：AP2 强（35×）+ 宽 RNA-binding 弱 + chromatin 弱阳；严格 PUF 未定（n=2）；curated GO 层因注释稀疏欠定。
- **CLM10-H-A（生命周期吸纳）假设**：富集模式确立，co-option 未确立（三极检验失败：配子体弱方向、zygote null、肝期检测偏倚排除）。
- **CLM10-H-C（ART 分区）假设/描述**：acute mRNA 正 vs chronic 蛋白负 vs 稳态转录 null 的分区描述；机制未立；无 translation×K13 公共数据（M4R-X01 缺口）。

### 3.3 已死亡（不许复活）

- 跨物种密码子应激程序（H3 跨物种版）；"Asn 是首个崩的供应"；"抗性正是通过这根轴实施"；"故意设计的瓶颈"；全基因组普遍密码子/氨基酸适应。

## 4. 现在还存在的漏洞（按严重度排序，每项附升级路线）

### 4.1 采样框偏差是最大漏洞（Figure 3 生死项）

**漏洞**：core-181 采样框排除 LCR 重/快进化蛋白（AP2/PUF 0 入选；core 成员中位长 259 vs 蛋白组其余 472 aa，Asn 0.056 vs 0.120）。在框内说"churn 不依赖 repeat"部分是被构造出来的。且项目约定的 LCR caller 在蛋白上 0%/0%（工具不可用）。
**已做**：D2b 全蛋白组归属已跑，SURVIVES 但只在框内。
**计算升级（必做）**：框外全蛋白组的 churn 归属（不限 core-OG）。若显示框外 churn 实为 LCR 富集，Figure 3 与 general-principle 强度必须降。——这是唯一能在投稿前用现有数据闭合的大漏洞。
**湿实验可闭合性**：低——这是序列层面问题，湿实验帮不上，只能靠更宽的采样框 + 更诚实的框外声明。

### 4.2 stage 混杂只有弱参考支撑（C 线措辞已收紧但仍受限）

**漏洞**：stage-adjusted 的 stage 参考 R² 仅 ~1e-4–3%；"stage 不是驱动者"只在该弱参考下成立，不能排除更强 stage 定义下的混杂。且 16/27 存活、6 组转弱，是如实计数。
**计算升级**：用每样本推断的发育年龄（inferred stage）替代名义时间重跑；需要高分辨 IDC 参考做 per-sample mapping。
**湿实验可闭合性**：中——如果做湿实验（§5），直接按真实 stage 同步收样（Percoll 同步化 + 时间点 + 形态镜检），把 stage 变成可测协变量，彻底绕过弱参考问题。

### 4.3 AP2 富集的方法边界（anchor 法 + IEA 循环）

**漏洞**：正式集实为 InterPro-anchored（anchor 发现用 product 文本），弱于人工 GO；真 GO 下全集阳性但 IEA 驱动，curated 非 IEA 层 n=10 欠定——是注释缺口，不是证伪。
**计算升级**：等 PlasmoDB 人工 GO 覆盖扩大后重跑非 IEA 层（T01 重激活条件：AP2 富集消失）。或者 orthology 级别构建严格 PUF 集（当前 n=2 不可判）。
**湿实验可闭合性**：低——注释问题，湿实验帮不上。

### 4.4 translation 层直接证据缺失（M4R-X01 缺口）

**漏洞**：acute DHA×K13 的 translation 层（Ribo/polysome/proteome/phospho/charging）无合格公共数据；scRNA 是转录层，不能代偿。
**计算升级**：无——公共数据就是没有。
**湿实验可闭合性**：高——这是§5 的核心靶点。做 §5.2 的 K13 配对 Ribo-seq/多聚核糖体谱，一次实验同时填这个缺口 + 区分 translation vs proteostasis。

### 4.5 SNP-level burden 未做（预算外）

**漏洞**：Pf8 SNP 水平负担检验需要 Zarr 流式基建，超 round-2 预算，已 deferred。
**计算升级**：建 Zarr 流式管线后补做（预算允许时）。
**湿实验可闭合性**：低。

### 4.6 阴性结果的统计效力边界

**漏洞**：D6 的 T2 AP2 是 4/4、p=0.0625（n=4，弱）；D3 的 MDR1/ABC 宽 domain 弱限制；C 线效应小（rho 0.06–0.23）。
**计算升级**：扩 AP2/PUF 宽 orthogroup 覆盖（目前 39 个宽 OG，AP2 只有 4 个 nsp≥8）；扩 transitions 数。
**湿实验可闭合性**：低-中。

## 5. 如果要补湿实验：最小验证怎么做（优先级 + 最小设计）

**原则**：只在以下两种情况下才值得开湿实验——(a) 一个 association 已经承重但缺因果闭合，且公共数据无法闭合（§4.4）；(b) 一个候选故事需要区分两个等价的竞争解释（当前：translation vs proteostasis vs generic stage effect）。**否则不开**（Gate D 原则：机制先于靶点，湿实验是条件性升级）。

### 5.1 第一优先：K13 配对的急性 translation 层（填 M4R-X01 缺口 + 区分竞争解释）

**要回答的问题**：急性 DHA 下，Asn 富集蛋白的输出变化，是**翻译层面的选择性维持**（translation），还是**蛋白稳定性/降解差异**（proteostasis），还是只是**发育阶段效应**？

**最小设计（一次实验同时区分）**：
- **模型**：等基因 K13 WT vs K13 C580Y（isogenic，消背景差异），DHA 急性 pulse（如 700nM×6h，沿用 RSA 时间窗）+ DMSO 对照；每组 ≥3 生物重复；同步化收样（Percoll + 镜检确认 stage，消 §4.2 漏洞）。
- **三层同时测**（同一份裂解物分用）：
  1. **Ribo-seq**（或简化 polysome profiling + 基因级 qPCR/测序定量）：测各基因翻译输出。
  2. **TMT 蛋白质组**：测蛋白稳态丰度。
  3. **靶向 tRNA 充电**（charge-seq 或 periodate-oxidation + northern/qPCR for Asn-tRNA(GUU)）：测 Asn 充电率（直接回应 L1c 的供应端）。
- **判读**：若 Ribo-seq 显示 Asn-top10 翻译输出维持而 TMT 未同步上升 → translation 选择性维持（支持 H-C 翻译臂）；若 TMT 上升但 Ribo-seq 无差异 → proteostasis；若两层的 stage 协变量吃掉全部信号 → generic stage effect（当前最诚实的默认）。**三种读法在实验开始前就写死**，不许事后调。
- **为什么是最小**：一个 pulse、两个基因型、三个层、≥3 重复；不跑全谱系、不跑长时间序列、不跑多药物。

### 5.2 第二优先：单基因座的 AP2 序列改造（验证"写在 Asn 上的调控"的因果性）

**要回答的问题**：AP2 蛋白的 Asn 富集，是不是功能必需的（不是装饰）？

**最小设计**：
- **靶向编辑**：选一个 AP2 候选（优先已有 perturbation 数据的），做两类编辑：①**同义重编码**（Asn 密码子 AAT→AAC，蛋白不变，密码子变）——区分密码子/翻译效应；②**保守替换**（Asn→Gln，IDR 性质保持，氨基酸变）——区分氨基酸化学效应。② 才是"written in Asn"的直接检验。
- **读出**：发育阶段转换率（gametocytogenesis 或 liver-stage development，取决于候选的功能域）+ parasite fitness。
- **对照**：未编辑 WT + 编辑无义对照（scramble 但保持 Asn 含量）——后者区分"是不是 Asn 本身"vs"是不是这个 IDR 的形状"。
- **预期与判读**：若 Asn→Gln 影响表型而同义重编码/scramble 无 → "written in Asn"有因果意义（化学臂）；若三者同效 → 是 IDR 结构/长度效应；若都无效 → 富集是装饰（当前最诚实的可能性之一）。

### 5.3 不做的湿实验（防漂移）

- 全谱系多基因型大panel——预算爆炸且不区分机制。
- 结构生物学/相分离成像——IDR≠phase separation，不许外推。
- 直接测"供应崩溃"动力学——L1c 已测过供应端，重复测不增加信息。

## 6. 创新点是什么（对着文献账逐项核）

### 6.1 我们的独特贡献（prior art 已逐项设防）

1. **branch-aware 的系统发育控制位点分析**：挂在进化分支上问"替换是否跟分支组成漂移同向"，LOO 7/7——这是对跨物种氨基酸差异的因果级控制，文献里没有同款。最接近的（Chaudhry 2018 / Battistuzzi 2016）只描述了 AT/LCR 耦合的存在，没有 phylogeny 控制的 churn + 架构分区 + 适应轴分离。
2. **组成轴 ≠ 适应轴的分离 + 群体层闭合**：14/14 耐药位点跨种保守零耦合 + Pf8 当代 34/34 分离——把"跨物种组成"和"当代适应"放进同一个框架并证明它们走不同的路。Hamilton 2017 MA 是单物种、无跨物种分解；Billows 2026 抬高了人群新颖性门槛，但我们的分离角度不是人群扫描。
3. **AP2 Asn 富集的双方法复现 + 诚实边界**：InterPro 与 regex 逐值一致（OR 35），同时自己量测了方法的 IEA/anchor 边界——候选级富集模式，不升机制。
4. **阴性结果的承重使用**：Gate C PIVOT、liver 检测偏倚排除、stage 弱参考声明、SNP deferred——不是删掉，而是作为叙事的刹车片和边界条件。

### 6.2 我们不做、也不能宣称的（防新颖性碰撞）

- **首次发现重复功能**：已被 2026 Nat Commun 抢先（NT-1）——我们 cite-and-complement。
- **动态 tRNA 修饰耐药**：归 Small-Saunders 2024 所有；Sinha-t6A 2025 我们 guard 为互补（不同位点 + 急性≠慢性分区）。
- **密码子适应**：M2 已降级。

## 7. 如果下一步只让我做一件事（主 Agent 建议）

按漏洞严重度 + 最小闭合成本排序：
1. **计算（投稿前必做）**：框外全蛋白组 churn 归属（§4.1）——一天工作量，决定 Figure 3 的强度上限。
2. **湿实验（若预算允许）**：§5.1 的 K13 配对 translation 层——一次实验同时填 M4R-X01 + 区分三个竞争解释 + 消 stage 混杂。
3. **如果都不做**：当前 D+A+C 按叙事版 HTML 打包投稿，T01 维持 superseded 但保留重激活条件（审稿要求家族级功能链 / AP2 在真 GO 下消失 / D2b 框外被证伪）。

## 8. 附：图版索引（与叙事版 HTML 共用 manuscript/figures/）

- Fig1a/b/c：组成地图与独立性；
- Fig2a/b/c：传导定量、branch-aware 桥、适应位点分离；
- Fig3a/b/c：架构分区（框内）、自然实验、Pf8 群体层；
- Fig4a/b/c：AP2 富集（含 GO 行）、反事实分层、单基因实证+阴性框；
- Fig5a/b + Fig5a_full_forest：ART 分区汇总与全森林（Extended）。

## 9. Claim 追溯（章节→claim→证据/决策）

- §1→PLAN.md §1–§2（未改写）
- §2.1→CLM01/02/03 + Gate A/B GO
- §2.2→Gate C PIVOT + D-016/D-017 + EVID-M3-001/003 + D-033
- §2.3→D-036 + L2 v2 + L1d/L1g/L1h + CLM08/09/10 注册
- §3→M5 终审 D-039 + 各线 ROUND 报告
- §4→对抗性审阅 §三 + L-006b/L-010 状态 + D2b verdict + stage/half-life 声明
- §5→M5 终审 §L6 + Gate D 原则 + L-007/L-009 最小验证思路
- §6→literature_anchor L026–L032 + M4R-L33–35 + NT-1/PA-1 设防
- 标题/结论→不超过 CANDIDATE/hypothesis



# Figures


## Fig1a

![](/home/huyudi/015_plasmo/manuscript/figures/Fig1a.png)


**Legend.** Figure 1 — Independent genome-composition transitions across Plasmodium (descriptive layer, CLM01/C1).  One-sentence claim: Plasmodium spans a 34.1pp genome-GC gradient (18.2–52.3%, n=18 species) that is not a single AT trajectory: ancestral reconstruction and leave-one-clade-out retain the gradient 7/7, with grade-A independence (low-GC Laverania+avian vs high-GC vivax+piroplasm poles).  (a) Genome GC (%) per species (n=18; frozen accessions per M1-01_composition_table.tsv; PlasmoDB-71 counterparts noted). Nuclear-chromosome GC tracks genome GC. (b) k-mer NJ species tree (M1-02_species_tree.nwk; branch lengths schematic — k-mer distances saturated) with tip genome-GC dots (blue→red) and squared-change-parsimony ancestral GC in red italics (17 internal nodes annotated from M1-02_ancestral_gc.tsv; root N0=43.6%). (c) Leave-one-clade-out: retained GC range per removed clade; gradient retained 7/7 (M1-02_leave_one_clade_out.tsv).  Stats units: species (n=18; genes do not multiply species N per I3). Sensitivity: compartment split (nuclear/mito/apicoplast in source table); N%/contig artifact audit in source README. Extended pointer: full QC tables in data/derived/WP1/M1-01_composition/ and M1-02_phylogeny/. No causal/adaptive language (C1 ceiling).


## Fig1b

![](/home/huyudi/015_plasmo/manuscript/figures/Fig1b.png)


**Legend.** Figure 1 — Independent genome-composition transitions across Plasmodium (descriptive layer, CLM01/C1).  One-sentence claim: Plasmodium spans a 34.1pp genome-GC gradient (18.2–52.3%, n=18 species) that is not a single AT trajectory: ancestral reconstruction and leave-one-clade-out retain the gradient 7/7, with grade-A independence (low-GC Laverania+avian vs high-GC vivax+piroplasm poles).  (a) Genome GC (%) per species (n=18; frozen accessions per M1-01_composition_table.tsv; PlasmoDB-71 counterparts noted). Nuclear-chromosome GC tracks genome GC. (b) k-mer NJ species tree (M1-02_species_tree.nwk; branch lengths schematic — k-mer distances saturated) with tip genome-GC dots (blue→red) and squared-change-parsimony ancestral GC in red italics (17 internal nodes annotated from M1-02_ancestral_gc.tsv; root N0=43.6%). (c) Leave-one-clade-out: retained GC range per removed clade; gradient retained 7/7 (M1-02_leave_one_clade_out.tsv).  Stats units: species (n=18; genes do not multiply species N per I3). Sensitivity: compartment split (nuclear/mito/apicoplast in source table); N%/contig artifact audit in source README. Extended pointer: full QC tables in data/derived/WP1/M1-01_composition/ and M1-02_phylogeny/. No causal/adaptive language (C1 ceiling).


## Fig1c

![](/home/huyudi/015_plasmo/manuscript/figures/Fig1c.png)


**Legend.** Figure 1 — Independent genome-composition transitions across Plasmodium (descriptive layer, CLM01/C1).  One-sentence claim: Plasmodium spans a 34.1pp genome-GC gradient (18.2–52.3%, n=18 species) that is not a single AT trajectory: ancestral reconstruction and leave-one-clade-out retain the gradient 7/7, with grade-A independence (low-GC Laverania+avian vs high-GC vivax+piroplasm poles).  (a) Genome GC (%) per species (n=18; frozen accessions per M1-01_composition_table.tsv; PlasmoDB-71 counterparts noted). Nuclear-chromosome GC tracks genome GC. (b) k-mer NJ species tree (M1-02_species_tree.nwk; branch lengths schematic — k-mer distances saturated) with tip genome-GC dots (blue→red) and squared-change-parsimony ancestral GC in red italics (17 internal nodes annotated from M1-02_ancestral_gc.tsv; root N0=43.6%). (c) Leave-one-clade-out: retained GC range per removed clade; gradient retained 7/7 (M1-02_leave_one_clade_out.tsv).  Stats units: species (n=18; genes do not multiply species N per I3). Sensitivity: compartment split (nuclear/mito/apicoplast in source table); N%/contig artifact audit in source README. Extended pointer: full QC tables in data/derived/WP1/M1-01_composition/ and M1-02_phylogeny/. No causal/adaptive language (C1 ceiling).


## Fig2a

![](/home/huyudi/015_plasmo/manuscript/figures/Fig2a.png)


**Legend.** Figure 2 — Base composition propagates quantitatively into codon and amino-acid space; composition-coupled turnover is branch-aware and distinct from resistance sites (CLM02/CLM03/CLM08, CANDIDATE/C2).  One-sentence claim: GC3 explains a median 84.5% of synonymous codon deviation and genome GC a median 78.7% of amino-acid CLR variation, yet composition-coupled sites follow branch GC/AT change far above background (fraction 0.791 vs 0.628, OR=2.25, Fisher p=9.5e-10, LOO OR range 1.71–4.35 significant 7/7), while 14/14 known resistance positions are composition-uncoupled at conserved sites.  (a) Background-explained variation (median across families; codon n=288 family×species rows, AA n=20; masked track; M2-01/M2-02 background_explained tables). Residuals are the minority. (b) Branch-aware concordance forest (n_events: coupled 316, uncoupled 8821, informative 13604; 95% CI; null 0.5 dashed). Symbols differ by stratum (circle/square/triangle) in addition to color. (c) Per-site separation: 14 known resistance positions (K13 539/580, CRT 76, DHFR 51/59/108/164, DHPS 436/437/540/581, MDR1 86/184/1246; OGs in table) all uncoupled (red x); matched-background coupling rises only at lower conservation (weak reading: distinct sequence spaces, not dual adaptive axes).  Stats units: parsimony sites (n=29,558; 646 ambiguous-root flagged) and branches; species are not independent replicates (branch-aware + LOO replace species-counting). Sensitivity: LCR masked/unmasked dual-track; uncertainty-site exclusion; terminal correlation retained as descriptive layer. Extended pointer: M4RD_D1_branch_substitutions.tsv (24,520 events), D1_LOO.tsv, D3 tables. Banned language respected: no dual-axes, no causal claims (association ceiling).


## Fig2b

![](/home/huyudi/015_plasmo/manuscript/figures/Fig2b.png)


**Legend.** Figure 2 — Base composition propagates quantitatively into codon and amino-acid space; composition-coupled turnover is branch-aware and distinct from resistance sites (CLM02/CLM03/CLM08, CANDIDATE/C2).  One-sentence claim: GC3 explains a median 84.5% of synonymous codon deviation and genome GC a median 78.7% of amino-acid CLR variation, yet composition-coupled sites follow branch GC/AT change far above background (fraction 0.791 vs 0.628, OR=2.25, Fisher p=9.5e-10, LOO OR range 1.71–4.35 significant 7/7), while 14/14 known resistance positions are composition-uncoupled at conserved sites.  (a) Background-explained variation (median across families; codon n=288 family×species rows, AA n=20; masked track; M2-01/M2-02 background_explained tables). Residuals are the minority. (b) Branch-aware concordance forest (n_events: coupled 316, uncoupled 8821, informative 13604; 95% CI; null 0.5 dashed). Symbols differ by stratum (circle/square/triangle) in addition to color. (c) Per-site separation: 14 known resistance positions (K13 539/580, CRT 76, DHFR 51/59/108/164, DHPS 436/437/540/581, MDR1 86/184/1246; OGs in table) all uncoupled (red x); matched-background coupling rises only at lower conservation (weak reading: distinct sequence spaces, not dual adaptive axes).  Stats units: parsimony sites (n=29,558; 646 ambiguous-root flagged) and branches; species are not independent replicates (branch-aware + LOO replace species-counting). Sensitivity: LCR masked/unmasked dual-track; uncertainty-site exclusion; terminal correlation retained as descriptive layer. Extended pointer: M4RD_D1_branch_substitutions.tsv (24,520 events), D1_LOO.tsv, D3 tables. Banned language respected: no dual-axes, no causal claims (association ceiling).


## Fig2c

![](/home/huyudi/015_plasmo/manuscript/figures/Fig2c.png)


**Legend.** Figure 2 — Base composition propagates quantitatively into codon and amino-acid space; composition-coupled turnover is branch-aware and distinct from resistance sites (CLM02/CLM03/CLM08, CANDIDATE/C2).  One-sentence claim: GC3 explains a median 84.5% of synonymous codon deviation and genome GC a median 78.7% of amino-acid CLR variation, yet composition-coupled sites follow branch GC/AT change far above background (fraction 0.791 vs 0.628, OR=2.25, Fisher p=9.5e-10, LOO OR range 1.71–4.35 significant 7/7), while 14/14 known resistance positions are composition-uncoupled at conserved sites.  (a) Background-explained variation (median across families; codon n=288 family×species rows, AA n=20; masked track; M2-01/M2-02 background_explained tables). Residuals are the minority. (b) Branch-aware concordance forest (n_events: coupled 316, uncoupled 8821, informative 13604; 95% CI; null 0.5 dashed). Symbols differ by stratum (circle/square/triangle) in addition to color. (c) Per-site separation: 14 known resistance positions (K13 539/580, CRT 76, DHFR 51/59/108/164, DHPS 436/437/540/581, MDR1 86/184/1246; OGs in table) all uncoupled (red x); matched-background coupling rises only at lower conservation (weak reading: distinct sequence spaces, not dual adaptive axes).  Stats units: parsimony sites (n=29,558; 646 ambiguous-root flagged) and branches; species are not independent replicates (branch-aware + LOO replace species-counting). Sensitivity: LCR masked/unmasked dual-track; uncertainty-site exclusion; terminal correlation retained as descriptive layer. Extended pointer: M4RD_D1_branch_substitutions.tsv (24,520 events), D1_LOO.tsv, D3 tables. Banned language respected: no dual-axes, no causal claims (association ceiling).


## Fig3a

![](/home/huyudi/015_plasmo/manuscript/figures/Fig3a.png)


**Legend.** Figure 3 — Composition-driven turnover partitions protein architecture within a quantified sampling frame (Route D; CLM08 CANDIDATE/C2, frame-internal qualifier mandatory).  One-sentence claim: Within conserved single-copy core genes, composition-coupled turnover falls in structured, domain-rich, Asn-depleted sequence with no LCR excess under length×OG-age matching (n_harbor=91 vs n_bg=5298; LCR MW p=0.10; Asn MW p=3.9e-33 depleted); two of three cross-lineage transitions are composition-consistent (T2+T3; T1 registered negative; CHROM null throughout); and 34/34 known resistance markers segregate in present-day Pf8 populations.  (a) Harbor-vs-matched medians (LCR fraction 0.0 vs 0.0; Asn 0.064 vs 0.119; polyN_max 2 vs 3; n_IPR 5 vs 1; MW asymptotic p annotated; bg hatched in addition to color). Right box: sampling-frame bias quantified (core-181 shorter/lower-Asn than rest of proteome) — "not repeat-confined" is frame-internal; D2b verdict SURVIVES. (b) Per-transition expected-sign fractions (AP2 circle, PUF_RNA square, CHROM triangle; n and sign-test p labeled; null 0.5 dashed). T2: AP2 4/4 (p=0.0625), PUF_RNA 21/27 (p=0.0030); T3: AP2 3/4, PUF_RNA 18/25 (p=0.022); T1: negative as registered (weak transition by construction: GC 0.19 vs 0.18). (c) Pf8 population layer (n≈24,409 samples; Zenodo 18681980): segregating vs fixed resistance markers; CNV amplification frequencies (GCH1/PM2-PM3 ~0.28, MDR1 0.20, CRT 0.086). SNP-level burden deferred (Zarr streaming beyond budget) — stated gap, not hidden.  Stats units: genes (D2b), per-transition sign tests (D6), population samples/markers (Pf8). Sensitivity: LCR caller dual-track (convention + 2 sensitivities); within-OG controls structurally impossible for fast-evolving drug genes (recorded limitation); terminal layer kept. Extended pointer: D2/D3/D4 tables, per-transition table (118 rows), M4RD_D2b_verdict.txt. Banned language respected: no repeat-first, no dual-axes, no mechanism (association ceiling); T1 negative shown, not rescued.


## Fig3b

![](/home/huyudi/015_plasmo/manuscript/figures/Fig3b.png)


**Legend.** Figure 3 — Composition-driven turnover partitions protein architecture within a quantified sampling frame (Route D; CLM08 CANDIDATE/C2, frame-internal qualifier mandatory).  One-sentence claim: Within conserved single-copy core genes, composition-coupled turnover falls in structured, domain-rich, Asn-depleted sequence with no LCR excess under length×OG-age matching (n_harbor=91 vs n_bg=5298; LCR MW p=0.10; Asn MW p=3.9e-33 depleted); two of three cross-lineage transitions are composition-consistent (T2+T3; T1 registered negative; CHROM null throughout); and 34/34 known resistance markers segregate in present-day Pf8 populations.  (a) Harbor-vs-matched medians (LCR fraction 0.0 vs 0.0; Asn 0.064 vs 0.119; polyN_max 2 vs 3; n_IPR 5 vs 1; MW asymptotic p annotated; bg hatched in addition to color). Right box: sampling-frame bias quantified (core-181 shorter/lower-Asn than rest of proteome) — "not repeat-confined" is frame-internal; D2b verdict SURVIVES. (b) Per-transition expected-sign fractions (AP2 circle, PUF_RNA square, CHROM triangle; n and sign-test p labeled; null 0.5 dashed). T2: AP2 4/4 (p=0.0625), PUF_RNA 21/27 (p=0.0030); T3: AP2 3/4, PUF_RNA 18/25 (p=0.022); T1: negative as registered (weak transition by construction: GC 0.19 vs 0.18). (c) Pf8 population layer (n≈24,409 samples; Zenodo 18681980): segregating vs fixed resistance markers; CNV amplification frequencies (GCH1/PM2-PM3 ~0.28, MDR1 0.20, CRT 0.086). SNP-level burden deferred (Zarr streaming beyond budget) — stated gap, not hidden.  Stats units: genes (D2b), per-transition sign tests (D6), population samples/markers (Pf8). Sensitivity: LCR caller dual-track (convention + 2 sensitivities); within-OG controls structurally impossible for fast-evolving drug genes (recorded limitation); terminal layer kept. Extended pointer: D2/D3/D4 tables, per-transition table (118 rows), M4RD_D2b_verdict.txt. Banned language respected: no repeat-first, no dual-axes, no mechanism (association ceiling); T1 negative shown, not rescued.


## Fig3c

![](/home/huyudi/015_plasmo/manuscript/figures/Fig3c.png)


**Legend.** Figure 3 — Composition-driven turnover partitions protein architecture within a quantified sampling frame (Route D; CLM08 CANDIDATE/C2, frame-internal qualifier mandatory).  One-sentence claim: Within conserved single-copy core genes, composition-coupled turnover falls in structured, domain-rich, Asn-depleted sequence with no LCR excess under length×OG-age matching (n_harbor=91 vs n_bg=5298; LCR MW p=0.10; Asn MW p=3.9e-33 depleted); two of three cross-lineage transitions are composition-consistent (T2+T3; T1 registered negative; CHROM null throughout); and 34/34 known resistance markers segregate in present-day Pf8 populations.  (a) Harbor-vs-matched medians (LCR fraction 0.0 vs 0.0; Asn 0.064 vs 0.119; polyN_max 2 vs 3; n_IPR 5 vs 1; MW asymptotic p annotated; bg hatched in addition to color). Right box: sampling-frame bias quantified (core-181 shorter/lower-Asn than rest of proteome) — "not repeat-confined" is frame-internal; D2b verdict SURVIVES. (b) Per-transition expected-sign fractions (AP2 circle, PUF_RNA square, CHROM triangle; n and sign-test p labeled; null 0.5 dashed). T2: AP2 4/4 (p=0.0625), PUF_RNA 21/27 (p=0.0030); T3: AP2 3/4, PUF_RNA 18/25 (p=0.022); T1: negative as registered (weak transition by construction: GC 0.19 vs 0.18). (c) Pf8 population layer (n≈24,409 samples; Zenodo 18681980): segregating vs fixed resistance markers; CNV amplification frequencies (GCH1/PM2-PM3 ~0.28, MDR1 0.20, CRT 0.086). SNP-level burden deferred (Zarr streaming beyond budget) — stated gap, not hidden.  Stats units: genes (D2b), per-transition sign tests (D6), population samples/markers (Pf8). Sensitivity: LCR caller dual-track (convention + 2 sensitivities); within-OG controls structurally impossible for fast-evolving drug genes (recorded limitation); terminal layer kept. Extended pointer: D2/D3/D4 tables, per-transition table (118 rows), M4RD_D2b_verdict.txt. Banned language respected: no repeat-first, no dual-axes, no mechanism (association ceiling); T1 negative shown, not rescued.


## Fig4a

![](/home/huyudi/015_plasmo/manuscript/figures/Fig4a.png)


**Legend.** Figure 4. Bounded regulatory-layer enrichment (Route A; CANDIDATE, enrichment != mechanism). (a) Asn top-decile membership odds ratios (95% CI) for regulator sets under InterPro-formal vs regex-baseline definitions. ApiAP2 replicates quantitatively across methods (OR 35.3 vs 34.3); chromatin flips to weak-positive only under the formal set; broad RNA-binding attenuates (3.08 -> 1.67) but stays significant; CCR4-NOT/transcription degenerate under formal sets (small-n, not refutation); sexual/gametocyte formal set is collinear with the AP2 set (do not double-count); strict PUF (n=2) is undetermined and is not evidence; PUF repeat-broad is null. L-010 true-GO check pending at build time. (b) Counterfactual stratification by LCR tertile: regulator Asn excess holds only in the high-LCR stratum (MW p=4.8e-07; zero/lowpos n.s.); length/LCR-adjusted logit OR=2.59 (1.69-3.97). The counterfactual is partial, not universal. (c) Left: single-locus perturbation proof (external positive control, GCN5 conditional knockdown source data): day-6 parasitaemia ~17x lower under +RAP; replicate line identity (LoxP/WT) not resolved from the snapshot layout, shown pooled by condition; see inventory. Middle: Pb ApiAP2 family map, descriptive only (11 KO-verified mutants, 9 transmission-blocked; small-n, not a test; Pb data not used as Pf expression evidence). Right: honest negative-report box for lifecycle state-consequence: gametocyte weak direction only (exploratory), zygote null, liver excluded as detection-driven. No shortlist; H-A stays hypothesis; family-wide lifecycle co-option is NOT claimed. N/units: universe n=5285 Pf genes, Asn top-decile q90=0.1837; Fisher OR + Woolf CI + BH q; MW asymptotic. Sensitivity/Extended: anchor-method caveat (formal sets are InterPro-anchored, weaker than curated GO) -> Extended Data + L-010.


## Fig4b

![](/home/huyudi/015_plasmo/manuscript/figures/Fig4b.png)


**Legend.** Figure 4. Bounded regulatory-layer enrichment (Route A; CANDIDATE, enrichment != mechanism). (a) Asn top-decile membership odds ratios (95% CI) for regulator sets under InterPro-formal vs regex-baseline definitions. ApiAP2 replicates quantitatively across methods (OR 35.3 vs 34.3); chromatin flips to weak-positive only under the formal set; broad RNA-binding attenuates (3.08 -> 1.67) but stays significant; CCR4-NOT/transcription degenerate under formal sets (small-n, not refutation); sexual/gametocyte formal set is collinear with the AP2 set (do not double-count); strict PUF (n=2) is undetermined and is not evidence; PUF repeat-broad is null. L-010 true-GO check pending at build time. (b) Counterfactual stratification by LCR tertile: regulator Asn excess holds only in the high-LCR stratum (MW p=4.8e-07; zero/lowpos n.s.); length/LCR-adjusted logit OR=2.59 (1.69-3.97). The counterfactual is partial, not universal. (c) Left: single-locus perturbation proof (external positive control, GCN5 conditional knockdown source data): day-6 parasitaemia ~17x lower under +RAP; replicate line identity (LoxP/WT) not resolved from the snapshot layout, shown pooled by condition; see inventory. Middle: Pb ApiAP2 family map, descriptive only (11 KO-verified mutants, 9 transmission-blocked; small-n, not a test; Pb data not used as Pf expression evidence). Right: honest negative-report box for lifecycle state-consequence: gametocyte weak direction only (exploratory), zygote null, liver excluded as detection-driven. No shortlist; H-A stays hypothesis; family-wide lifecycle co-option is NOT claimed. N/units: universe n=5285 Pf genes, Asn top-decile q90=0.1837; Fisher OR + Woolf CI + BH q; MW asymptotic. Sensitivity/Extended: anchor-method caveat (formal sets are InterPro-anchored, weaker than curated GO) -> Extended Data + L-010.


## Fig4c

![](/home/huyudi/015_plasmo/manuscript/figures/Fig4c.png)


**Legend.** Figure 4. Bounded regulatory-layer enrichment (Route A; CANDIDATE, enrichment != mechanism). (a) Asn top-decile membership odds ratios (95% CI) for regulator sets under InterPro-formal vs regex-baseline definitions. ApiAP2 replicates quantitatively across methods (OR 35.3 vs 34.3); chromatin flips to weak-positive only under the formal set; broad RNA-binding attenuates (3.08 -> 1.67) but stays significant; CCR4-NOT/transcription degenerate under formal sets (small-n, not refutation); sexual/gametocyte formal set is collinear with the AP2 set (do not double-count); strict PUF (n=2) is undetermined and is not evidence; PUF repeat-broad is null. L-010 true-GO check pending at build time. (b) Counterfactual stratification by LCR tertile: regulator Asn excess holds only in the high-LCR stratum (MW p=4.8e-07; zero/lowpos n.s.); length/LCR-adjusted logit OR=2.59 (1.69-3.97). The counterfactual is partial, not universal. (c) Left: single-locus perturbation proof (external positive control, GCN5 conditional knockdown source data): day-6 parasitaemia ~17x lower under +RAP; replicate line identity (LoxP/WT) not resolved from the snapshot layout, shown pooled by condition; see inventory. Middle: Pb ApiAP2 family map, descriptive only (11 KO-verified mutants, 9 transmission-blocked; small-n, not a test; Pb data not used as Pf expression evidence). Right: honest negative-report box for lifecycle state-consequence: gametocyte weak direction only (exploratory), zygote null, liver excluded as detection-driven. No shortlist; H-A stays hypothesis; family-wide lifecycle co-option is NOT claimed. N/units: universe n=5285 Pf genes, Asn top-decile q90=0.1837; Fisher OR + Woolf CI + BH q; MW asymptotic. Sensitivity/Extended: anchor-method caveat (formal sets are InterPro-anchored, weaker than curated GO) -> Extended Data + L-010.


## Fig5a

![](/home/huyudi/015_plasmo/manuscript/figures/Fig5a.png)


**Legend.** Figure 5. ART-stress partition, second context (Route C; CONDITIONAL, not MOA). (a) Per-contrast Spearman rho between gene Asn fraction and response (95% Fisher-z CI; * q<0.05), acute (blue circles) vs chronic (red squares) arms. GSE151189 rows shown stage-adjusted (STAGEADJ preferred; 16/27 raw-significant survive adjusted and positive). Stage reference is weak (R2~1e-4-3%), so "stage is not the driver" holds only under this reference; length attenuation is large but nonzero; acute evidence is mRNA-layer only (microarray + scRNA), cross-dataset independence downgraded half a grade. Chronic protein (Mok2021) runs opposite-sign; chronic transcript (GSE59099) is null; genotype interaction (DiD) is weak. Mobilization/conservation language is hypothesis only. (b) Left: dTE interaction upgrade under rule 2 — direction agrees with ratio-based dTE (sign agreement 0.924, spearman 0.996) but only 1 gene at q<0.05, so WEAK SUPPORT retained and ratio-based dTE stays DISCOVERY/SUPPORTING. Right: standing M4R-X01 gap — no suitable public acute DHA x K13 translation-layer data exist; scRNA is transcript-only and does not close it. C stays a partition figure, not a mechanism figure. N/units: per-contrast gene n in panel table; q from validation table (BH per analysis). Sensitivity/Extended: no-Schizont tracks and full 142-contrast table in Extended Data.


## Fig5b

![](/home/huyudi/015_plasmo/manuscript/figures/Fig5b.png)


**Legend.** Figure 5. ART-stress partition, second context (Route C; CONDITIONAL, not MOA). (a) Per-contrast Spearman rho between gene Asn fraction and response (95% Fisher-z CI; * q<0.05), acute (blue circles) vs chronic (red squares) arms. GSE151189 rows shown stage-adjusted (STAGEADJ preferred; 16/27 raw-significant survive adjusted and positive). Stage reference is weak (R2~1e-4-3%), so "stage is not the driver" holds only under this reference; length attenuation is large but nonzero; acute evidence is mRNA-layer only (microarray + scRNA), cross-dataset independence downgraded half a grade. Chronic protein (Mok2021) runs opposite-sign; chronic transcript (GSE59099) is null; genotype interaction (DiD) is weak. Mobilization/conservation language is hypothesis only. (b) Left: dTE interaction upgrade under rule 2 — direction agrees with ratio-based dTE (sign agreement 0.924, spearman 0.996) but only 1 gene at q<0.05, so WEAK SUPPORT retained and ratio-based dTE stays DISCOVERY/SUPPORTING. Right: standing M4R-X01 gap — no suitable public acute DHA x K13 translation-layer data exist; scRNA is transcript-only and does not close it. C stays a partition figure, not a mechanism figure. N/units: per-contrast gene n in panel table; q from validation table (BH per analysis). Sensitivity/Extended: no-Schizont tracks and full 142-contrast table in Extended Data.


## Fig5a_full_forest

![](/home/huyudi/015_plasmo/manuscript/figures/Fig5a_full_forest.png)


**Legend.** Figure 5. ART-stress partition, second context (Route C; CONDITIONAL, not MOA). (a) Per-contrast Spearman rho between gene Asn fraction and response (95% Fisher-z CI; * q<0.05), acute (blue circles) vs chronic (red squares) arms. GSE151189 rows shown stage-adjusted (STAGEADJ preferred; 16/27 raw-significant survive adjusted and positive). Stage reference is weak (R2~1e-4-3%), so "stage is not the driver" holds only under this reference; length attenuation is large but nonzero; acute evidence is mRNA-layer only (microarray + scRNA), cross-dataset independence downgraded half a grade. Chronic protein (Mok2021) runs opposite-sign; chronic transcript (GSE59099) is null; genotype interaction (DiD) is weak. Mobilization/conservation language is hypothesis only. (b) Left: dTE interaction upgrade under rule 2 — direction agrees with ratio-based dTE (sign agreement 0.924, spearman 0.996) but only 1 gene at q<0.05, so WEAK SUPPORT retained and ratio-based dTE stays DISCOVERY/SUPPORTING. Right: standing M4R-X01 gap — no suitable public acute DHA x K13 translation-layer data exist; scRNA is transcript-only and does not close it. C stays a partition figure, not a mechanism figure. N/units: per-contrast gene n in panel table; q from validation table (BH per analysis). Sensitivity/Extended: no-Schizont tracks and full 142-contrast table in Extended Data.



# Supplementary Tables


## TableS1_D1_concordance.tsv — D1 branch-aware concordance

| stratum | n_events | n_concordant | frac | ci_lo | ci_hi | binom_p |
| --- | --- | --- | --- | --- | --- | --- |
| all_branches | 24520 | 13704 | 0.5588907014681892 | 0.5526670814883219 | 0.5650958712784256 | 4.48146742049739e-76 |
| informative_branches | 13604 | 7927 | 0.5826962658041752 | 0.5743875782710643 | 0.5909582617334224 | 3.220405927382791e-83 |
| informative_coupled | 316 | 250 | 0.7911392405063291 | 0.7429557863027209 | 0.832328979656429 | 2.369582470187817e-26 |
| informative_uncoupled | 8821 | 5538 | 0.627819975059517 | 0.6176786962919285 | 0.637849969514379 | 1.1368041495472738e-128 |
| informative_1nt | 9472 | 5173 | 0.5461359797297297 | 0.5360928107545629 | 0.5561417407389951 | 2.8000013758277044e-19 |
| informative_gen2nt | 4132 | 2754 | 0.6665053242981607 | 0.6519810363069953 | 0.6807202934537963 | 1.9654827686019305e-103 |


## TableS1b_D1_LOO.tsv — D1 branch leave-one-out

| left_out | n_species | n_events | n_sites | threshold | all_branches_n | all_branches_frac | all_branches_lo | all_branches_hi | all_branches_p | informative_branches_n | informative_branches_frac | informative_branches_lo | informative_branches_hi | informative_branches_p | informative_coupled_n | informative_coupled_frac | informative_coupled_lo | informative_coupled_hi | informative_coupled_p | informative_uncoupled_n | informative_uncoupled_frac | informative_uncoupled_lo | informative_uncoupled_hi | informative_uncoupled_p | Fisher_OR | Fisher_p |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| vivax | 12 | 20617 | 26959 | 0.03717258767772147 | 20617 | 0.5579861279526604 | 0.5511968334282025 | 0.5647538171989767 | 2.4060260955370184e-62 | 11235 | 0.5939474855362706 | 0.5848358477085831 | 0.6029948981220193 | 1.0956856892723872e-88 | 136 | 0.7573529411764706 | 0.6788806041686491 | 0.8216857509073797 | 1.419571074706475e-09 | 7288 | 0.6458562019758507 | 0.6348017993510691 | 0.6567569201910559 | 9.705445203194895e-139 | 1.711461330964199 | 0.006478225016424621 |
| malariae | 15 | 23372 | 29137 | 0.0375638450847271 | 23372 | 0.5612271093616293 | 0.554855508162647 | 0.5675785863599649 | 2.4454324045753632e-78 | 13550 | 0.5852398523985239 | 0.5769211553020118 | 0.5935102299961608 | 5.555061749729923e-88 | 314 | 0.802547770700637 | 0.7549744976082314 | 0.8428075422506166 | 2.5430603200591756e-28 | 8807 | 0.6298399000794822 | 0.6197008845358318 | 0.6398656930769128 | 1.509820413808049e-132 | 2.3887367190634867 | 7.610262042747705e-11 |
| rodent | 13 | 21688 | 28155 | 0.038046367448600515 | 21688 | 0.5589265953522685 | 0.5523085978637486 | 0.565523721177342 | 1.40771160349661e-67 | 11268 | 0.590610578629748 | 0.5815019086911047 | 0.59965748588378 | 7.84502848207333e-83 | 241 | 0.8132780082987552 | 0.7593033056704304 | 0.8574219567032838 | 1.1806312481144184e-23 | 7535 | 0.643662906436629 | 0.6327785100875676 | 0.6544008889330989 | 4.279386227489725e-139 | 2.41127147766323 | 1.9735915643733904e-08 |
| avian | 15 | 23322 | 29154 | 0.03882316275354791 | 23322 | 0.5779092702169625 | 0.5715581647041057 | 0.5842347135258484 | 1.3555399453666122e-125 | 11761 | 0.5912762520193862 | 0.5823631208131211 | 0.600129773948845 | 1.2494874441030014e-87 | 305 | 0.8262295081967214 | 0.7797175347762325 | 0.8646257143197442 | 3.464961081146188e-32 | 7822 | 0.6343646126310406 | 0.6236279667419323 | 0.6449693429710575 | 2.9149794108495876e-126 | 2.7405261116561337 | 5.589131924575424e-13 |
| piroplasm | 13 | 12097 | 20754 | 0.14213315012463626 | 12097 | 0.5955195503017277 | 0.5867444705864909 | 0.6042339816890929 | 1.5930429556947634e-98 | 6335 | 0.6640883977900552 | 0.6523612899723508 | 0.6756166169057766 | 4.969933150589339e-153 | 96 | 0.9270833333333334 | 0.8570680669280426 | 0.9642328737620105 | 3.259152789937673e-19 | 4463 | 0.7452386287250728 | 0.7322477595117639 | 0.7578076747351098 | 1.2423356316811973e-245 | 4.34640494802852 | 1.2171866570798713e-05 |
| coccidian | 15 | 20002 | 27312 | 0.03881693482754056 | 20002 | 0.5148485151484852 | 0.5079200845223744 | 0.5217712432346143 | 2.7501847338269136e-05 | 10656 | 0.5212087087087087 | 0.5117177589450315 | 0.5306843720584912 | 1.2447879689799748e-05 | 240 | 0.7875 | 0.7314251575172397 | 0.8345160108477973 | 7.938753608792669e-20 | 6798 | 0.5641365107384525 | 0.5523157585660838 | 0.5758848158195355 | 3.753622128407106e-26 | 2.863241046092492 | 1.386638613531959e-12 |
| laverania | 13 | 22196 | 28242 | 0.037666445097389256 | 22196 | 0.579744098035682 | 0.5732371223440392 | 0.58622347488833 | 2.960774776334958e-125 | 13145 | 0.5943704830734119 | 0.5859500886242076 | 0.6027357344676151 | 2.2577564185624028e-104 | 257 | 0.7937743190661478 | 0.7401567550751222 | 0.8387386358441175 | 4.517076045182682e-22 | 8573 | 0.6402659512422723 | 0.6300459803880502 | 0.6503602707895395 | 1.3261463293767433e-150 | 2.1625962044156926 | 1.9480495696387733e-07 |


## TableS2_D2b_matched.tsv — D2b harbor vs matched background

| feature | n_harbor | n_bg | median_harbor | median_bg | MW_p | median_strat_diff | n_strata |
| --- | --- | --- | --- | --- | --- | --- | --- |
| lcr_frac | 91 | 5298 | 0.0 | 0.0 | 0.1036393355383494 | 0.0 | 8 |
| asn_frac | 91 | 5298 | 0.0642201834862385 | 0.1193120193255768 | 3.88286261226388e-33 | -0.031035275520439738 | 8 |
| polyN_max | 91 | 5298 | 2.0 | 3.0 | 2.381782968042223e-18 | -0.375 | 8 |
| n_ipr | 91 | 5298 | 5.0 | 1.0 | 5.678961033122527e-28 | 3.0 | 8 |


## TableS2b_D2b_framebias.tsv — D2b sampling-frame bias quant

| feature | median_core | median_noncore | MW_p |
| --- | --- | --- | --- |
| lcr_frac | 0.0 | 0.0 | 0.02211004327361117 |
| asn_frac | 0.0560747663551401 | 0.12012289069924065 | 4.865026083405028e-75 |
| length | 259.0 | 472.0 | 7.763755354216337e-23 |


## TableS3_D6formal_summary.tsv — D6 formal summary

| transition | set | n_genes | median_diff_Asn | n_expected_sign | signtest_p | expected_sign |
| --- | --- | --- | --- | --- | --- | --- |
| T1_AT_Lav | AP2 | 4 | 0.0003354978354978312 | 2 | 0.6875 | 1 |
| T1_AT_Lav | CHROM | 5 | 0.0 | 0 | 1.0 | 1 |
| T1_AT_Lav | PUF_RNA | 25 | -0.00023310023310023353 | 9 | 0.9461239278316498 | 1 |
| T2_GC_vivax | AP2 | 4 | -0.01879494993578278 | 4 | 0.0625 | -1 |
| T2_GC_vivax | CHROM | 6 | 0.0 | 2 | 0.890625 | -1 |
| T2_GC_vivax | PUF_RNA | 27 | -0.016355140186915886 | 21 | 0.0029623061418533325 | -1 |
| T3_rodent | AP2 | 4 | 0.018295983455309917 | 3 | 0.3125 | 1 |
| T3_rodent | CHROM | 6 | 0.0 | 2 | 0.890625 | 1 |
| T3_rodent | PUF_RNA | 25 | 0.007741771865846087 | 18 | 0.021642625331878662 | 1 |


## TableS3b_D6formal_pertransition.tsv — D6 per-transition tests

| gene | set | og | transition | n_A | n_B | mean_Asn_A | mean_Asn_B | diff_Asn | diff_LCR | status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| PF3D7_0934400 | AP2 | OG_004127 | T1_AT_Lav | 2 | 1 | 0.055 | 0.05 | 0.0049999999999999975 | 0.0 | OK |
| PF3D7_0934400 | AP2 | OG_004127 | T2_GC_vivax | 4 | 4 | 0.04125443937087772 | 0.05694538478898669 | -0.01569094541810897 | 0.0 | OK |
| PF3D7_0934400 | AP2 | OG_004127 | T3_rodent | 3 | 1 | 0.05529225908372828 | 0.06190476190476191 | -0.006612502821033629 | 0.0 | OK |
| PF3D7_0611200 | AP2 | OG_005156 | T1_AT_Lav | 2 | 1 | 0.1079136690647482 | 0.10071942446043165 | 0.007194244604316544 | 0.0 | OK |
| PF3D7_0611200 | AP2 | OG_005156 | T2_GC_vivax | 4 | 4 | 0.09069139099104853 | 0.10627620097014831 | -0.015584809979099778 | 0.0 | OK |
| PF3D7_0611200 | AP2 | OG_005156 | T3_rodent | 3 | 1 | 0.1087216248506571 | 0.0989399293286219 | 0.0097816955220352 | 0.0 | OK |
| PF3D7_1115500 | AP2 | OG_000484 | T1_AT_Lav | 2 | 1 | 0.10822510822510822 | 0.11255411255411256 | -0.004329004329004335 | 0.0 | OK |
| PF3D7_1115500 | AP2 | OG_000484 | T2_GC_vivax | 4 | 4 | 0.07305957678547273 | 0.11219427161553315 | -0.03913469483006042 | 0.0 | OK |
| PF3D7_1115500 | AP2 | OG_000484 | T3_rodent | 3 | 1 | 0.12005649717514126 | 0.08860759493670886 | 0.0314489022384324 | 0.0 | OK |
| PF3D7_1305200 | AP2 | OG_002139 | T1_AT_Lav | 2 | 1 | 0.11128048780487805 | 0.11585365853658537 | -0.004573170731707321 | 0.0 | OK |
| PF3D7_1305200 | AP2 | OG_002139 | T2_GC_vivax | 4 | 4 | 0.08857019487111442 | 0.11046914932457101 | -0.02189895445345659 | 0.0 | OK |
| PF3D7_1305200 | AP2 | OG_002139 | T3_rodent | 3 | 1 | 0.11717171717171716 | 0.09036144578313253 | 0.026810271388584633 | 0.0 | OK |
| PF3D7_0823200 | PUF_RNA | OG_001763 | T1_AT_Lav | 2 | 1 | 0.06643356643356643 | 0.06666666666666667 | -0.00023310023310023353 | 0.0 | OK |
| PF3D7_0823200 | PUF_RNA | OG_001763 | T2_GC_vivax | 4 | 4 | 0.057974147871408144 | 0.07640229515229516 | -0.018428147280887014 | 0.0 | OK |
| PF3D7_0823200 | PUF_RNA | OG_001763 | T3_rodent | 3 | 1 | 0.0787037037037037 | 0.0694980694980695 | 0.009205634205634208 | 0.0 | OK |
| PF3D7_1406000 | PUF_RNA | OG_000838 | T1_AT_Lav | 2 | 1 | 0.09345794392523364 | 0.09345794392523364 | 0.0 | 0.0 | OK |
| PF3D7_1406000 | PUF_RNA | OG_000838 | T2_GC_vivax | 4 | 4 | 0.0834665262922547 | 0.08644859813084112 | -0.0029820718385864114 | 0.0 | OK |
| PF3D7_1406000 | PUF_RNA | OG_000838 | T3_rodent | 3 | 1 | 0.08722741433021806 | 0.08411214953271028 | 0.0031152647975077885 | 0.0 | OK |
| PF3D7_1347500 | PUF_RNA | OG_002473 | T1_AT_Lav | 2 | 1 | 0.0913978494623656 | 0.08870967741935484 | 0.0026881720430107503 | 0.0 | OK |
| PF3D7_1347500 | PUF_RNA | OG_002473 | T2_GC_vivax | 4 | 4 | 0.08176943699731903 | 0.07625697122621898 | 0.005512465771100053 | 0.0 | OK |
| PF3D7_1347500 | PUF_RNA | OG_002473 | T3_rodent | 3 | 1 | 0.07486631016042782 | 0.08042895442359249 | -0.005562644263164676 | 0.0 | OK |
| PF3D7_0617300 | PUF_RNA | OG_005213 | T1_AT_Lav | 2 | 1 | 0.08992805755395683 | 0.09352517985611511 | -0.0035971223021582788 | 0.0 | OK |
| PF3D7_0617300 | PUF_RNA | OG_005213 | T2_GC_vivax | 4 | 4 | 0.09892086330935251 | 0.09172661870503597 | 0.007194244604316544 | 0.0 | OK |
| PF3D7_0617300 | PUF_RNA | OG_005213 | T3_rodent | 3 | 1 | 0.09352517985611512 | 0.08633093525179857 | 0.0071942446043165575 | 0.0 | OK |
| PF3D7_1346300 | PUF_RNA | OG_002463 | T1_AT_Lav | 2 | 1 | 0.08095421664890359 | 0.08095238095238096 | 1.8356965226368205e-06 | 0.0 | OK |
| PF3D7_1346300 | PUF_RNA | OG_002463 | T2_GC_vivax | 4 | 4 | 0.058960990582982986 | 0.05953909728789154 | -0.000578106704908557 | 0.0 | OK |
| PF3D7_1346300 | PUF_RNA | OG_002463 | T3_rodent | 3 | 1 | 0.05499521914808303 | 0.07317073170731707 | -0.01817551255923404 | 0.0 | OK |
| PF3D7_1024200 | PUF_RNA | OG_000215 | T1_AT_Lav | 2 | 1 | 0.084375 | 0.0875 | -0.003124999999999989 | 0.0 | OK |
| PF3D7_1024200 | PUF_RNA | OG_000215 | T2_GC_vivax | 4 | 4 | 0.054041701104711654 | 0.09762986256696948 | -0.04358816146225783 | 0.0 | OK |
| PF3D7_1024200 | PUF_RNA | OG_000215 | T3_rodent | 3 | 1 | 0.10082304526748971 | 0.0880503144654088 | 0.012772730802080912 | 0.0 | OK |
| PF3D7_1360100 | PUF_RNA | OG_002573 | T1_AT_Lav | 2 | 1 | 0.12777777777777777 | 0.12568306010928962 | 0.0020947176684881497 | 0.0 | OK |
| PF3D7_1360100 | PUF_RNA | OG_002573 | T2_GC_vivax | 4 | 4 | 0.12569060773480661 | 0.1111111111111111 | 0.01457949662369551 | 0.0 | OK |
| PF3D7_1360100 | PUF_RNA | OG_002573 | T3_rodent | 3 | 1 | 0.11111111111111112 | 0.1111111111111111 | 1.3877787807814457e-17 | 0.0 | OK |
| PF3D7_0814200 | PUF_RNA | OG_001820 | T1_AT_Lav | 2 | 1 | 0.028225806451612902 | 0.028225806451612902 | 0.0 | 0.0 | OK |
| PF3D7_0814200 | PUF_RNA | OG_001820 | T2_GC_vivax | 3 | 4 | 0.026886261693472534 | 0.031122784595587098 | -0.004236522902114563 | 0.0 | OK |
| PF3D7_0814200 | PUF_RNA | OG_001820 | T3_rodent | 3 | 1 | 0.03305822756204862 | 0.02531645569620253 | 0.007741771865846087 | 0.0 | OK |
| PF3D7_0728900 | PUF_RNA | OG_004418 | T1_AT_Lav | 2 | 1 | 0.12041016926562709 | 0.12448979591836734 | -0.004079626652740254 | 0.0 | OK |
| PF3D7_0728900 | PUF_RNA | OG_004418 | T2_GC_vivax | 4 | 4 | 0.07673343066039975 | 0.12648108438895017 | -0.04974765372855042 | 0.028368794326241134 | OK |
| PF3D7_0728900 | PUF_RNA | OG_004418 | T3_rodent | 3 | 1 | 0.13626049347098118 | 0.09714285714285714 | 0.03911763632812404 | 0.0 | OK |
| PF3D7_0205700 | PUF_RNA | OG_001957 | T1_AT_Lav | 2 | 1 | 0.14143877802414387 | 0.14634146341463414 | -0.004902685390490269 | 0.0 | OK |
| PF3D7_0205700 | PUF_RNA | OG_001957 | T2_GC_vivax | 4 | 4 | 0.11680297940703958 | 0.12807713247595337 | -0.011274153068913786 | 0.0 | OK |
| PF3D7_0205700 | PUF_RNA | OG_001957 | T3_rodent | 3 | 1 | 0.12686707094354757 | 0.13170731707317074 | -0.004840246129623166 | 0.0 | OK |
| PF3D7_1135500 | PUF_RNA | OG_004963 | T1_AT_Lav | 2 | 1 | 0.1055486539247988 | 0.10555555555555556 | -6.901630756755983e-06 | 0.0 | OK |
| PF3D7_1135500 | PUF_RNA | OG_004963 | T2_GC_vivax | 4 | 4 | 0.08260869565217391 | 0.11227770177838578 | -0.029669006126211872 | 0.0 | OK |
| PF3D7_1135500 | PUF_RNA | OG_004963 | T3_rodent | 3 | 1 | 0.10852713178294575 | 0.12352941176470589 | -0.015002279981760139 | 0.0 | OK |
| PF3D7_1467500 | PUF_RNA | OG_001390 | T1_AT_Lav | 2 | 1 | 0.09389140271493213 | 0.09502262443438914 | -0.0011312217194570096 | 0.0 | OK |
| PF3D7_1467500 | PUF_RNA | OG_001390 | T2_GC_vivax | 4 | 4 | 0.052598493910272506 | 0.07332743879202572 | -0.02072894488175321 | 0.0 | OK |
| PF3D7_1467500 | PUF_RNA | OG_001390 | T3_rodent | 3 | 1 | 0.074391336187256 | 0.07013574660633484 | 0.0042555895809211625 | 0.0 | OK |
| PF3D7_1415300 | PUF_RNA | OG_000927 | T1_AT_Lav | 2 | 1 | 0.10682492581602374 | 0.10979228486646884 | -0.002967359050445109 | 0.0 | OK |
| PF3D7_1415300 | PUF_RNA | OG_000927 | T2_GC_vivax | 4 | 4 | 0.09412811336416613 | 0.10796319313808692 | -0.013835079773920789 | 0.0 | OK |
| PF3D7_1415300 | PUF_RNA | OG_000927 | T3_rodent | 3 | 1 | 0.10823663846983018 | 0.10714285714285714 | 0.0010937813269730479 | 0.0 | OK |
| PF3D7_0417100 | PUF_RNA | OG_003564 | T1_AT_Lav | 2 | 1 | 0.10797665369649806 | 0.11673151750972763 | -0.008754863813229569 | 0.0 | OK |
| PF3D7_0417100 | PUF_RNA | OG_003564 | T2_GC_vivax | 3 | 4 | 0.060951867134837634 | 0.09719938720148097 | -0.036247520066643335 | 0.0 | OK |
| PF3D7_0417100 | PUF_RNA | OG_003564 | T3_rodent | 3 | 1 | 0.10544459356332729 | 0.07246376811594203 | 0.03298082544738526 | 0.0 | OK |
| PF3D7_1207500 | PUF_RNA | OG_002707 | T1_AT_Lav | 2 | 1 | 0.08241758241758242 | 0.07692307692307693 | 0.005494505494505489 | 0.0 | OK |
| PF3D7_1207500 | PUF_RNA | OG_002707 | T2_GC_vivax | 3 | 4 | 0.07875457875457875 | 0.07768569674647022 | 0.0010688820081085348 | 0.0 | OK |
| PF3D7_1207500 | PUF_RNA | OG_002707 | T3_rodent | 3 | 1 | 0.08148148148148147 | 0.06629834254143646 | 0.015183138940045013 | 0.0 | OK |
| PF3D7_0812500 | PUF_RNA | OG_001833 | T1_AT_Lav | 2 | 1 | 0.1495822446308916 | 0.14600840336134455 | 0.0035738412695470456 | 0.0 | OK |
| PF3D7_0812500 | PUF_RNA | OG_001833 | T2_GC_vivax | 4 | 3 | 0.09003224784576097 | 0.12659431030513082 | -0.036562062459369854 | 0.06330067822155237 | OK |
| PF3D7_0812500 | PUF_RNA | OG_001833 | T3_rodent | 2 | 1 | 0.13438524876142627 | 0.11101243339253997 | 0.0233728153688863 | 0.0 | OK |
| PF3D7_1454000 | PUF_RNA | OG_001267 | T1_AT_Lav | 2 | 1 | 0.15166340508806261 | 0.16237623762376238 | -0.010712832535699768 | 0.0 | OK |
| PF3D7_1454000 | PUF_RNA | OG_001267 | T2_GC_vivax | 4 | 3 | 0.10791475995773936 | 0.14114455800418754 | -0.03322979804644818 | 0.0 | OK |
| PF3D7_1454000 | PUF_RNA | OG_001267 | T3_rodent | 2 | 1 | 0.14818548387096775 | 0.12706270627062707 | 0.021122777600340675 | 0.0 | OK |
| PF3D7_1330800 | PUF_RNA | OG_002347 | T1_AT_Lav | 2 | 1 | 0.16571318832376486 | 0.16220735785953178 | 0.003505830464233084 | 0.04480651731160898 | OK |
| PF3D7_1330800 | PUF_RNA | OG_002347 | T2_GC_vivax | 3 | 3 | 0.14706896874609385 | 0.14631451843779725 | 0.0007544503082966003 | -0.09775967413441955 | OK |
| PF3D7_1330800 | PUF_RNA | OG_002347 | T3_rodent | 2 | 1 | 0.13959466552151922 | 0.1597542242703533 | -0.020159558748834072 | 0.08961303462321792 | OK |
| PF3D7_1317300 | PUF_RNA | OG_002244 | T1_AT_Lav | 2 | 1 | 0.10623691873691873 | 0.11148648648648649 | -0.005249567749567752 | 0.0 | OK |
| PF3D7_1317300 | PUF_RNA | OG_002244 | T2_GC_vivax | 3 | 4 | 0.05951038851083674 | 0.08357829710726242 | -0.02406790859642568 | 0.0 | OK |
| PF3D7_1317300 | PUF_RNA | OG_002244 | T3_rodent | 3 | 1 | 0.08812770616632659 | 0.06993006993006994 | 0.018197636236256656 | 0.0 | OK |
| PF3D7_1320900 | PUF_RNA | OG_004612 | T1_AT_Lav | 2 | 0 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_1320900 | PUF_RNA | OG_004612 | T2_GC_vivax | 4 | 4 | 0.10910556543335309 | 0.0886376722440945 | 0.020467893189258596 | 0.0 | OK |
| PF3D7_1320900 | PUF_RNA | OG_004612 | T3_rodent | 3 | 1 | 0.08618356299212598 | 0.096 | -0.009816437007874018 | 0.0 | OK |
| PF3D7_1006200 | PUF_RNA | OG_000053 | T1_AT_Lav | 1 | 1 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_1006200 | PUF_RNA | OG_000053 | T2_GC_vivax | 4 | 4 | 0.028037383177570093 | 0.04439252336448598 | -0.016355140186915886 | 0.0 | OK |
| PF3D7_1006200 | PUF_RNA | OG_000053 | T3_rodent | 3 | 1 | 0.04672897196261682 | 0.037383177570093455 | 0.009345794392523366 | 0.0 | OK |
| PF3D7_0929200 | PUF_RNA | OG_004559 | T1_AT_Lav | 1 | 1 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0929200 | PUF_RNA | OG_004559 | T2_GC_vivax | 4 | 4 | 0.14491486455892735 | 0.16904871235602498 | -0.024133847797097624 | 0.0 | OK |
| PF3D7_0929200 | PUF_RNA | OG_004559 | T3_rodent | 3 | 1 | 0.16845878136200718 | 0.1708185053380783 | -0.002359723976071121 | 0.0 | OK |
| PF3D7_0319500 | PUF_RNA | OG_003394 | T1_AT_Lav | 2 | 1 | 0.06273086855611128 | 0.0627062706270627 | 2.4597929048575384e-05 | 0.0 | OK |
| PF3D7_0319500 | PUF_RNA | OG_003394 | T2_GC_vivax | 3 | 2 | 0.06094975711903868 | 0.07167577413479054 | -0.010726017015751856 | 0.0 | OK |
| PF3D7_0319500 | PUF_RNA | OG_003394 | T3_rodent | 1 | 1 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_1326300 | PUF_RNA | OG_002318 | T1_AT_Lav | 2 | 1 | 0.13008624832080587 | 0.1189083820662768 | 0.011177866254529073 | 0.11098779134295227 | OK |
| PF3D7_1326300 | PUF_RNA | OG_002318 | T2_GC_vivax | 2 | 4 | 0.053045391384900414 | 0.10948615654498008 | -0.05644076516007966 | -0.06659267480577137 | OK |
| PF3D7_1326300 | PUF_RNA | OG_002318 | T3_rodent | 3 | 1 | 0.11175034116210587 | 0.1026936026936027 | 0.009056738468503175 | 0.08879023307436183 | OK |
| PF3D7_0916700 | PUF_RNA | OG_004005 | T1_AT_Lav | 2 | 0 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0916700 | PUF_RNA | OG_004005 | T2_GC_vivax | 3 | 4 | 0.10392156862745099 | 0.12577027089403453 | -0.021848702266583545 | 0.0 | OK |
| PF3D7_0916700 | PUF_RNA | OG_004005 | T3_rodent | 3 | 1 | 0.13054510014786933 | 0.11144578313253012 | 0.019099317015339207 | 0.0 | OK |
| PF3D7_1360900 | PUF_RNA | OG_004819 | T1_AT_Lav | 2 | 1 | 0.13624678663239073 | 0.13881748071979436 | -0.0025706940874036244 | 0.0 | OK |
| PF3D7_1360900 | PUF_RNA | OG_004819 | T2_GC_vivax | 4 | 1 | 0.1283213757012435 | 0.13554987212276215 | -0.007228496421518654 | 0.0 | OK |
| PF3D7_1360900 | PUF_RNA | OG_004819 | T3_rodent | 0 | 1 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0811900 | PUF_RNA | OG_004454 | T1_AT_Lav | 2 | 1 | 0.12018244226392569 | 0.12521008403361344 | -0.005027641769687757 | 0.0 | OK |
| PF3D7_0811900 | PUF_RNA | OG_004454 | T2_GC_vivax | 1 | 4 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0811900 | PUF_RNA | OG_004454 | T3_rodent | 3 | 1 | 0.13736049189225583 | 0.10490111779879621 | 0.032459374093459614 | 0.02037178507766743 | OK |
| PF3D7_0606100 | PUF_RNA | OG_005110 | T1_AT_Lav | 2 | 1 | 0.13773804714987384 | 0.1343804537521815 | 0.0033575933976923433 | 0.0 | OK |
| PF3D7_0606100 | PUF_RNA | OG_005110 | T2_GC_vivax | 3 | 2 | 0.0800412353551935 | 0.1332990178347848 | -0.05325778247959129 | 0.08852005532503458 | OK |
| PF3D7_0606100 | PUF_RNA | OG_005110 | T3_rodent | 1 | 1 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0820100 | PUF_RNA | OG_004475 | T1_AT_Lav | 2 | 1 | 0.09179433611884866 | 0.09705882352941177 | -0.0052644874105631095 | 0.0 | OK |
| PF3D7_0820100 | PUF_RNA | OG_004475 | T2_GC_vivax | 1 | 3 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0820100 | PUF_RNA | OG_004475 | T3_rodent | 3 | 0 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_0320900 | CHROM | OG_003406 | T1_AT_Lav | 2 | 1 | 0.0189873417721519 | 0.0189873417721519 | 0.0 | 0.0 | OK |
| PF3D7_0320900 | CHROM | OG_003406 | T2_GC_vivax | 4 | 4 | 0.0189873417721519 | 0.0189873417721519 | 0.0 | 0.0 | OK |
| PF3D7_0320900 | CHROM | OG_003406 | T3_rodent | 3 | 1 | 0.0189873417721519 | 0.0189873417721519 | 0.0 | 0.0 | OK |
| PF3D7_0617800 | CHROM | OG_005217 | T1_AT_Lav | 2 | 1 | 0.06060606060606061 | 0.06060606060606061 | 0.0 | 0.0 | OK |
| PF3D7_0617800 | CHROM | OG_005217 | T2_GC_vivax | 4 | 4 | 0.06203007518796992 | 0.06203007518796992 | 0.0 | 0.0 | OK |
| PF3D7_0617800 | CHROM | OG_005217 | T3_rodent | 3 | 1 | 0.06015037593984962 | 0.06766917293233082 | -0.007518796992481203 | 0.0 | OK |
| PF3D7_1105100 | CHROM | OG_000397 | T1_AT_Lav | 2 | 1 | 0.017094017094017096 | 0.017094017094017096 | 0.0 | 0.0 | OK |
| PF3D7_1105100 | CHROM | OG_000397 | T2_GC_vivax | 4 | 4 | 0.017057800956106043 | 0.023305084745762712 | -0.006247283789656669 | 0.0 | OK |
| PF3D7_1105100 | CHROM | OG_000397 | T3_rodent | 3 | 1 | 0.025423728813559324 | 0.01694915254237288 | 0.008474576271186442 | 0.0 | OK |
| PF3D7_0714000 | CHROM | OG_001577 | T1_AT_Lav | 2 | 1 | 0.032520325203252036 | 0.032520325203252036 | 0.0 | 0.0 | OK |
| PF3D7_0714000 | CHROM | OG_001577 | T2_GC_vivax | 4 | 4 | 0.03265470671235638 | 0.032520325203252036 | 0.00013438150910434488 | 0.0 | OK |
| PF3D7_0714000 | CHROM | OG_001577 | T3_rodent | 3 | 1 | 0.032520325203252036 | 0.032520325203252036 | 0.0 | 0.0 | OK |
| PF3D7_0925700 | CHROM | OG_004073 | T1_AT_Lav | 2 | 1 | 0.051224944320712694 | 0.051224944320712694 | 0.0 | 0.0 | OK |
| PF3D7_0925700 | CHROM | OG_004073 | T2_GC_vivax | 4 | 4 | 0.04930066593145571 | 0.04977628635346756 | -0.0004756204220118451 | 0.0 | OK |
| PF3D7_0925700 | CHROM | OG_004073 | T3_rodent | 3 | 1 | 0.04996271439224459 | 0.049217002237136466 | 0.0007457121551081233 | 0.0 | OK |
| PF3D7_1105000 | CHROM | OG_000396 | T1_AT_Lav | 2 | 0 |  |  |  |  | UNRESOLVED_coverage |
| PF3D7_1105000 | CHROM | OG_000396 | T2_GC_vivax | 4 | 4 | 0.019417475728155338 | 0.019417475728155338 | 0.0 | 0.0 | OK |
| PF3D7_1105000 | CHROM | OG_000396 | T3_rodent | 3 | 1 | 0.019417475728155338 | 0.019417475728155338 | 0.0 | 0.0 | OK |


## TableS3c_D3_multidim.tsv — D3 multidim matched controls

| gene | pf_pos | scope | n_bg | frac_bg_coupled | note |
| --- | --- | --- | --- | --- | --- |
| PF3D7_1343700 | 580 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_1343700 | 580 | domain_shared | 531 | 0.0449438202247191 |  |
| PF3D7_1343700 | 580 | known_gene_ess_bin | 1 |  | MIS=0.122 |
| PF3D7_1343700 | 539 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_1343700 | 539 | domain_shared | 531 | 0.0449438202247191 |  |
| PF3D7_1343700 | 539 | known_gene_ess_bin | 1 |  | MIS=0.122 |
| PF3D7_0709000 | 76 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0709000 | 76 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0709000 | 76 | known_gene_ess_bin | 1 |  | MIS=0.127 |
| PF3D7_0417200 | 51 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 51 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 51 | known_gene_ess_bin | 1 |  | MIS=0.125 |
| PF3D7_0417200 | 59 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 59 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 59 | known_gene_ess_bin | 1 |  | MIS=0.125 |
| PF3D7_0417200 | 108 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 108 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 108 | known_gene_ess_bin | 1 |  | MIS=0.125 |
| PF3D7_0417200 | 164 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 164 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0417200 | 164 | known_gene_ess_bin | 1 |  | MIS=0.125 |
| PF3D7_0810800 | 436 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 436 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 436 | known_gene_ess_bin | 1 |  | MIS=0.195 |
| PF3D7_0810800 | 437 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 437 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 437 | known_gene_ess_bin | 1 |  | MIS=0.195 |
| PF3D7_0810800 | 540 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 540 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 540 | known_gene_ess_bin | 1 |  | MIS=0.195 |
| PF3D7_0810800 | 581 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 581 | domain_shared | 0 |  | n<5:uninformative |
| PF3D7_0810800 | 581 | known_gene_ess_bin | 1 |  | MIS=0.195 |
| PF3D7_0523000 | 86 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0523000 | 86 | domain_shared | 5056 | 0.020496894409937887 |  |
| PF3D7_0523000 | 86 | known_gene_ess_bin | 1 |  | MIS=0.128 |
| PF3D7_0523000 | 184 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0523000 | 184 | domain_shared | 5056 | 0.020496894409937887 |  |
| PF3D7_0523000 | 184 | known_gene_ess_bin | 1 |  | MIS=0.128 |
| PF3D7_0523000 | 1246 | within_OG | 0 |  | n<5:uninformative |
| PF3D7_0523000 | 1246 | domain_shared | 5056 | 0.020496894409937887 |  |
| PF3D7_0523000 | 1246 | known_gene_ess_bin | 1 |  | MIS=0.128 |


## TableS3d_Pf8boundary.tsv — Pf8 population boundary

| marker | gene | pos | ref | call_rate | n_alleles | segregating | res_freq_all | res_freq_strict | n_strict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| crt_72[C] | crt | 72 | C | 0.9986070711622762 | 3 | True | 0.01723076923076923 | 0.02289002557544757 | 15640 |
| crt_74[M] | crt | 74 | M | 0.9978696382481872 | 2 | True | 0.5353286529539762 | 0.6056707629288275 | 15624 |
| crt_75[N] | crt | 75 | N | 0.9978696382481872 | 6 | True | 0.5399679763517674 | 0.612224 | 15625 |
| crt_76[K] | crt | 76 | K | 0.9980335122290959 | 4 | True | 0.5579820204425106 | 0.636736 | 15625 |
| crt_93[T] | crt | 93 | T | 0.9988528821336392 | 4 | True | 0.04425577293794348 | 0.055299833780846436 | 15642 |
| crt_97[H] | crt | 97 | H | 0.9987299766479577 | 4 | True | 0.017228648781688408 | 0.02142080695696656 | 15639 |
| crt_218[I] | crt | 218 | I | 0.9710352738743906 | 2 | True | 0.02594717745337946 | 0.035062281684571275 | 15173 |
| crt_220[A] | crt | 220 | A | 0.9546069072882953 | 4 | True | 0.5443114029440796 | 0.6257878503419606 | 14914 |
| crt_271[Q] | crt | 271 | Q | 0.9922569544020649 | 3 | True | 0.5338976052848885 | 0.6055512622359609 | 15528 |
| crt_326[N] | crt | 326 | N | 0.9847597197754926 | 3 | True | 0.26958439073095647 | 0.3491258854877494 | 15387 |
| crt_333[T] | crt | 333 | T | 0.9932401982875169 | 6 | True | 0.0532090414123082 | 0.06553897607409313 | 15548 |
| crt_353[G] | crt | 353 | G | 0.998811913638412 | 3 | True | 0.004142739950779327 | 0.004857782038990093 | 15645 |
| crt_356[I] | crt | 356 | I | 0.9959850874677373 | 3 | True | 0.38932993295216156 | 0.4606258817493908 | 15594 |
| crt_371[R] | crt | 371 | R | 0.9964357409152362 | 3 | True | 0.5085108132554889 | 0.5639300134589502 | 15603 |
| dhfr_16[N] | dhfr | 16 | N | 0.9957392764963743 | 2 | True | 1.0 | 1.0 | 15591 |
| dhfr_51[N] | dhfr | 51 | N | 0.9878323569175305 | 2 | True | 0.8329047777040478 | 0.8064432656229784 | 15458 |
| dhfr_59[C] | dhfr | 59 | C | 0.986767176041624 | 3 | True | 0.9155526031719671 | 0.9097510373443983 | 15424 |
| dhfr_108[S] | dhfr | 108 | S | 0.9830390429759515 | 3 | True | 0.9493644509272765 | 0.9485978267941961 | 15369 |
| dhfr_164[I] | dhfr | 164 | I | 0.9889794747838911 | 3 | True | 0.20853355426677714 | 0.26270857586340707 | 15462 |
| dhfr_306[S] | dhfr | 306 | S | 0.9893072227457086 | 2 | True | 0.011015405002484678 | 0.012926577042399173 | 15472 |
| dhps_436[S] | dhps | 436 | S | 0.9863984595845795 | 13 | True | 0.32923536985504837 | 0.2762262683859263 | 15433 |
| dhps_437[G] | dhps | 437 | G | 0.988938506288664 | 2 | True | 0.20170678155681676 | 0.17067494181536075 | 15468 |
| dhps_540[K] | dhps | 540 | K | 0.9851284362325372 | 7 | True | 0.46560758546119935 | 0.5010074748131297 | 15385 |
| dhps_581[A] | dhps | 581 | A | 0.9850464992420829 | 2 | True | 0.22246714357012146 | 0.26767709620763674 | 15373 |
| dhps_613[A] | dhps | 613 | A | 0.987504608955713 | 4 | True | 0.05642217059409227 | 0.03787043641787173 | 15421 |
| exo_415[E] | exo | 415 | E | 0.9986480396575034 | 3 | True | 0.10895963242533639 | 0.1415872609835646 | 15637 |
| mdr1_86[N] | mdr1 | 86 | N | 0.9863165225941252 | 5 | True | 0.12191069574247144 | 0.09129222589638851 | 15423 |
| mdr1_184[Y] | mdr1 | 184 | Y | 0.969191691589168 | 2 | True | 0.5550154288371306 | 0.5047247736734289 | 15133 |
| mdr1_1034[S] | mdr1 | 1034 | S | 0.949977467327625 | 3 | True | 0.002803174055545972 | 0.0039824502193722576 | 14815 |
| mdr1_1042[N] | mdr1 | 1042 | N | 0.9558769306403376 | 2 | True | 0.016800960054860276 | 0.02265567397278638 | 14919 |
| mdr1_1226[F] | mdr1 | 1226 | F | 0.9694784710557581 | 2 | True | 0.04377958079783638 | 0.05232251253628926 | 15156 |
| mdr1_1246[D] | mdr1 | 1246 | D | 0.968249416198943 | 2 | True | 0.03905390539053905 | 0.028826446280991735 | 15125 |
| fd_193[D] | fd | 193 | D | 0.9972141423245524 | 4 | True | 0.22201224271804773 | 0.2879011080509832 | 15613 |
| mdr2_484[T] | mdr2 | 484 | T | 0.9926666393543365 | 2 | True | 0.2306231943871234 | 0.2948635427394439 | 15536 |


## TableS4_interpro_enrichment.tsv — InterPro-formal vs regex enrichment

| set | n_set | n_in_top10 | OR | lo | hi | p | q_bh |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ApiAP2 [InterPro-formal] | 24 | 19 | 35.32289628180039 | 13.133548926714361 | 95.00151167799918 | 2.0612137326919728e-15 | 1.2367282396151836e-14 |
| ApiAP2 [regex-baseline] | 28 | 22 | 34.27755905511811 | 13.834373623306735 | 84.92983396065532 | 1.5059735283918084e-17 | 2.710752351105255e-16 |
| PUF [InterPro-formal] | 0 | 0 |  |  |  | 1.0 | 1.0 |
| PUF [regex-baseline] | 2 | 1 | 8.986767485822305 | 0.5612657933338439 | 143.89259207214033 | 0.19052787480027528 | 0.263807826646535 |
| RNA_binding_broad [InterPro-formal] | 467 | 70 | 1.6704632570364693 | 1.273105222264953 | 2.1918435682358335 | 0.00036590116192321607 | 0.00094088870208827 |
| RNA_binding_broad [regex-baseline] | 201 | 49 | 3.084951854688697 | 2.2058950845412872 | 4.314315767980534 | 1.848189145097964e-09 | 8.316851152940838e-09 |
| chromatin_reg [InterPro-formal] | 176 | 31 | 1.975122659111326 | 1.3257088632385647 | 2.942659302288289 | 0.0019070844467618709 | 0.0038141688935237417 |
| chromatin_reg [regex-baseline] | 131 | 19 | 1.5413929270338271 | 0.9397689681357475 | 2.5281662153870106 | 0.10279804758744368 | 0.18503648565739864 |
| CCR4_NOT [InterPro-formal] | 8 | 2 | 2.9981060606060606 | 0.6035801309206648 | 14.892206502775467 | 0.18765619611882595 | 0.263807826646535 |
| CCR4_NOT [regex-baseline] | 12 | 6 | 9.062977099236642 | 2.9124658055452515 | 28.202066353843634 | 0.0005387093320433309 | 0.0012120959970974947 |
| transcription_reg [InterPro-formal] | 4 | 0 | 0.0 | 0.0 | inf | 1.0 | 1.0 |
| transcription_reg [regex-baseline] | 89 | 26 | 3.8420256991685564 | 2.4108561681878244 | 6.122788106503758 | 2.832114475706655e-07 | 1.0195612112543957e-06 |
| sexual_gametocyte [InterPro-formal] | 24 | 19 | 35.32289628180039 | 13.133548926714361 | 95.00151167799918 | 2.0612137326919728e-15 | 1.2367282396151836e-14 |
| sexual_gametocyte [regex-baseline] | 45 | 13 | 3.711254835589942 | 1.9355074935730472 | 7.116176248567911 | 0.00031864143476774843 | 0.00094088870208827 |
| proteostasis [InterPro-formal] | 318 | 35 | 1.1173216261555483 | 0.7769953862465696 | 1.6067117493522833 | 0.5629580289930994 | 0.6755496347917194 |
| proteostasis [regex-baseline] | 166 | 19 | 1.1655417548624147 | 0.7165510246939774 | 1.8958699876368696 | 0.5117915262107234 | 0.6580176765566444 |
| PUF_strict [InterPro-formal] | 2 | 1 | 8.986767485822305 | 0.5612657933338439 | 143.89259207214033 | 0.19052787480027528 | 0.263807826646535 |
| PUF_repeatbroad [InterPro-formal] | 79 | 9 | 1.156155744447491 | 0.5741636080136768 | 2.3280752850972948 | 0.7041603746781677 | 0.7921804215129387 |


## TableS5_lifecycle_tripole.tsv — Lifecycle tripole

| pole | contrast | test | n | stat | p | med_hit | med_rest | q_bh |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| gametocyte | female-vs-male | spearman(asn_frac,logFC) | 4949 | -0.08811359983505142 | 5.324520675117036e-10 |  |  | 2.7383249186316184e-09 |
| gametocyte | female-vs-male | MW Asn-top10 vs rest | 4949 | 1104239.0 | 0.0367840071302104 | -0.133840545817722 | 0.12400936664398345 | 0.08107137874062806 |
| gametocyte | female-vs-male | MW AP2(IPR001471) vs rest | 4949 | 46246.0 | 0.06565455045512707 | -0.9541167772421546 | 0.0958545487985481 | 0.13130910091025413 |
| gametocyte | female-vs-male | MW PUF-strict vs rest | 4949 |  |  | -0.3936593641223042 | 0.0938753453279199 | 1.0 |
| zygote | 12_h-vs-0_h | spearman(asn_frac,logFC) | 4953 | 0.022823225821892296 | 0.10826368142319932 |  |  | 0.19487462656175877 |
| zygote | 12_h-vs-0_h | MW Asn-top10 vs rest | 4953 | 1185924.5 | 0.37400100319396645 | -0.5778252860118425 | -0.6115704131459836 | 0.5385614445993118 |
| zygote | 12_h-vs-0_h | MW AP2(IPR001471) vs rest | 4953 | 66464.5 | 0.2948800758156398 | -0.27743936527126045 | -0.6060203721398612 | 0.46155142301578406 |
| zygote | 12_h-vs-0_h | MW PUF-strict vs rest | 4953 |  |  | 1.1192176936339533 | -0.6051121750024206 | 1.0 |
| zygote | 20_h-vs-0_h | spearman(asn_frac,logFC) | 4953 | -0.013755086656838035 | 0.33311889411372875 |  |  | 0.4996783411705931 |
| zygote | 20_h-vs-0_h | MW Asn-top10 vs rest | 4953 | 1113304.5 | 0.14422065324746958 | -1.4806178681463162 | -1.198591453933481 | 0.24723540556709073 |
| zygote | 20_h-vs-0_h | MW AP2(IPR001471) vs rest | 4953 | 64523.0 | 0.4416348524633441 | -0.7457674596931789 | -1.2313750934107575 | 0.6114944111030918 |
| zygote | 20_h-vs-0_h | MW PUF-strict vs rest | 4953 |  |  | 1.8545580865819726 | -1.2313750934107575 | 1.0 |
| zygote | 20_h-vs-12_h | spearman(asn_frac,logFC) | 4953 | -0.025740773292688295 | 0.07007655574593846 |  |  | 0.13277663193967287 |
| zygote | 20_h-vs-12_h | MW Asn-top10 vs rest | 4953 | 1090982.0 | 0.028212668655423732 | -0.0638326104224145 | 0.0 | 0.07254686225680387 |
| zygote | 20_h-vs-12_h | MW AP2(IPR001471) vs rest | 4953 | 57388.5 | 0.8001106747999132 | 0.0 | 0.0 | 1.0 |
| zygote | 20_h-vs-12_h | MW PUF-strict vs rest | 4953 |  |  | 0.7353403929480193 | 0.0 | 1.0 |
| liver | GFP-vs-NoGFP(all-days) | spearman(asn_frac,logFC) | 5285 | 0.1034058519915711 | 4.852412790752976e-14 |  |  | 4.3671715116776784e-13 |
| liver | GFP-vs-NoGFP(all-days) | MW Asn-top10 vs rest | 5285 | 1479890.0 | 4.15753673517232e-11 | 4.340512175782142 | 3.3181945289800567 | 2.9934264493240706e-10 |
| liver | GFP-vs-NoGFP(all-days) | MW AP2(IPR001471) vs rest | 5285 | 79885.5 | 0.024659186689623133 | 5.134024525416986 | 3.427844529817188 | 0.06828697852511022 |
| liver | GFP-vs-NoGFP(all-days) | MW PUF-strict vs rest | 5285 |  |  | 4.538886644223367 | 3.4410363267179704 | 1.0 |
| liver | GFP-vs-NoGFP(Day2) | spearman(asn_frac,logFC) | 5285 | -0.028499729160386223 | 0.03828370662751881 |  |  | 0.08107137874062806 |
| liver | GFP-vs-NoGFP(Day2) | MW Asn-top10 vs rest | 5285 | 1335535.5 | 0.003144986499997385 | 0.0 | 0.0 | 0.011321951399990586 |
| liver | GFP-vs-NoGFP(Day2) | MW AP2(IPR001471) vs rest | 5285 | 76596.5 | 0.01856799344351542 | 0.0 | 0.0 | 0.055703980330546264 |
| liver | GFP-vs-NoGFP(Day2) | MW PUF-strict vs rest | 5285 |  |  | 3.7913755718440836 | 0.0 | 1.0 |
| liver | GFP-vs-NoGFP(Day4) | spearman(asn_frac,logFC) | 5285 | 0.06466447679252033 | 2.541263371387083e-06 |  |  | 1.0165053485548332e-05 |
| liver | GFP-vs-NoGFP(Day4) | MW Asn-top10 vs rest | 5285 | 1460122.5 | 3.776612253369067e-10 | 4.308366938811153 | 2.663754484512854 | 2.26596735202144e-09 |
| liver | GFP-vs-NoGFP(Day4) | MW AP2(IPR001471) vs rest | 5285 | 78581.5 | 0.030707329880435945 | 6.731694140712138 | 2.663754484512854 | 0.07369759171304627 |
| liver | GFP-vs-NoGFP(Day4) | MW PUF-strict vs rest | 5285 |  |  | 4.481592401248015 | 2.663754484512854 | 1.0 |
| liver | GFP-vs-NoGFP(Day6) | spearman(asn_frac,logFC) | 5285 | 0.10453175524001666 | 2.5697711928827693e-14 |  |  | 3.083725431459323e-13 |
| liver | GFP-vs-NoGFP(Day6) | MW Asn-top10 vs rest | 5285 | 1430996.0 | 2.88696720964258e-07 | 7.455842184525327 | 6.781407398671224 | 1.299135244339161e-06 |
| liver | GFP-vs-NoGFP(Day6) | MW AP2(IPR001471) vs rest | 5285 | 71570.0 | 0.25782995245197526 | 8.022450963968595 | 6.836274176594336 | 0.42190355855777767 |
| liver | GFP-vs-NoGFP(Day6) | MW PUF-strict vs rest | 5285 |  |  | 5.343691959578002 | 6.8388687667356525 | 1.0 |
| liver | GFP-vs-Day0naive | spearman(asn_frac,logFC) | 5285 | 0.11357522090245917 | 1.2158516075768224e-16 |  |  | 2.1885328936382805e-15 |
| liver | GFP-vs-Day0naive | MW Asn-top10 vs rest | 5285 | 1538224.0 | 6.859784123477446e-17 | 4.793115949462549 | 3.54134556639976 | 2.1885328936382805e-15 |
| liver | GFP-vs-Day0naive | MW AP2(IPR001471) vs rest | 5285 | 82271.0 | 0.010269607411065789 | 5.75685909558502 | 3.6454487341801673 | 0.033609624254397126 |
| liver | GFP-vs-Day0naive | MW PUF-strict vs rest | 5285 |  |  | 4.538886644223367 | 3.653267903522335 | 1.0 |


## TableS6_C_partition_asnfrac.tsv — C acute-vs-chronic partition

| contrast | feature | n | spearman_rho | spearman_p | lenadj_rho | lenadj_p | MW_AsnTop10_p | median_Top10 | median_rest | arm |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GSE151189_DHA-CTRL_Cam3II_R539T_8.0h | asn_frac | 5148 | 0.061921442002344146 | 8.747407815762815e-06 | 0.09445512058436326 | 1.1162633284618187e-11 | 0.4481614227192363 | 0.0044720365 | 0.0352948699999999 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_8.0h_noSchiz | asn_frac | 3626 | 0.0827129229646601 | 6.116730429729019e-07 | 0.10360925292971844 | 4.009726140821932e-10 | 0.41321870183243836 | 0.0021069857499999 | 0.0356181946666667 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_8.0h_STAGEADJ | asn_frac | 4767 | 0.043634317621753196 | 0.0025841521490319683 | 0.07896535024597685 | 4.784534353969062e-08 | 0.26887084303942177 | -0.055108747113928325 | -0.017064190302403567 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_8.0h_STAGEADJ_noSchiz | asn_frac | 3245 | 0.05754154840623016 | 0.0010406386674144295 | 0.08202254694031828 | 2.889942681472108e-06 | 0.13884094823965568 | -0.06154890874881944 | -0.016646027477871316 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_16.0h | asn_frac | 5092 | 0.07989475492286142 | 1.137236066748159e-08 | 0.06250558992982161 | 8.060927324494109e-06 | 0.6211684078221343 | 0.00904963191666655 | 0.00207768408333335 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_16.0h_noSchiz | asn_frac | 3586 | 0.10041270244480822 | 1.6788551898708674e-09 | 0.06752015514625115 | 5.198182021978382e-05 | 0.3236961056192875 | -0.00430204 | -0.013283012 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_16.0h_STAGEADJ | asn_frac | 4721 | 0.06532391579890463 | 7.05622751041226e-06 | 0.05431370212876887 | 0.00018873455146145232 | 0.8699541783557101 | -0.030148304749099554 | -0.03784285977598861 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_16.0h_STAGEADJ_noSchiz | asn_frac | 3215 | 0.09081139854336039 | 2.499328393774904e-07 | 0.05841733653302187 | 0.000920170149618051 | 0.5953986382797799 | -0.0507997151705149 | -0.05115705173968275 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_24.0h | asn_frac | 5142 | 0.06183189131672449 | 9.122430263333295e-06 | 0.0330619928078954 | 0.017746151298754318 | 0.18217954069624254 | 0.0356345166666666 | 0.0226336245 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_24.0h_noSchiz | asn_frac | 3622 | 0.06558797327303954 | 7.80938427202513e-05 | 0.03360797453506542 | 0.04312434113281282 | 0.9184920392092093 | -0.0563839396666663 | -0.03669409 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_24.0h_STAGEADJ | asn_frac | 4763 | 0.018298481796446192 | 0.20672050850843873 | 0.01776289683610942 | 0.22032204946265527 | 0.7956950242622156 | -0.01527016938947361 | -0.022811467220013502 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_24.0h_STAGEADJ_noSchiz | asn_frac | 3243 | 0.04132667879980907 | 0.01859547000114161 | 0.013595894221411072 | 0.43893827035238087 | 0.8327382860250888 | -0.07827901952346236 | -0.06748069032330267 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_32.0h | asn_frac | 5126 | 0.03539114882666656 | 0.011275467495223275 | 0.020019984941926448 | 0.15181592352505757 | 0.848927258800227 | 0.0059514867499999 | -0.00396679 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_32.0h_noSchiz | asn_frac | 3615 | 0.04512131124210146 | 0.006660537723268702 | 0.004459100197521959 | 0.7886905186794408 | 0.9722832415709143 | -0.0006174129722222501 | -0.0025666435 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_32.0h_STAGEADJ | asn_frac | 4743 | 0.06574997854716122 | 5.846245438447601e-06 | 0.0170719243340635 | 0.23979087502915464 | 0.06768900032960891 | 0.024094002209463228 | -0.029894136614617213 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_32.0h_STAGEADJ_noSchiz | asn_frac | 3232 | 0.0591858730181507 | 0.0007615429916214038 | 0.003950188363175411 | 0.8223800099038565 | 0.9447317007082461 | -0.04969618377585483 | -0.050309901653793676 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_40.0h | asn_frac | 5123 | -0.09993475342129134 | 7.551998213389653e-13 | -0.051119109410410056 | 0.00025186363044798674 | 4.070002960206714e-08 | -0.0086668321111111 | 0.07830946853703699 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_40.0h_noSchiz | asn_frac | 3612 | -0.08717462232436968 | 1.5428937189331354e-07 | -0.04487331219945476 | 0.006990267931178026 | 7.417851012058124e-06 | 0.0798177595833333 | 0.16827246216666658 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_40.0h_STAGEADJ | asn_frac | 4743 | -0.09566456764700056 | 4.062665990317142e-11 | -0.07688003076807232 | 1.1508427692542861e-07 | 1.0319955291766212e-05 | -0.028613568885701005 | 0.0422739471073614 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_40.0h_STAGEADJ_noSchiz | asn_frac | 3232 | -0.13647235139940855 | 6.5943133938777695e-15 | -0.08776301485495154 | 5.819844101110302e-07 | 7.548165027299205e-05 | 0.04660302782235618 | 0.11649386817132441 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_48.0h | asn_frac | 5108 | -0.008886628415287439 | 0.5254369862436061 | -0.0566800094340033 | 5.0521168520116725e-05 | 0.9983852893507226 | -0.0144671914999999 | -0.0264663851666666 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_48.0h_noSchiz | asn_frac | 3602 | -0.009379536946720436 | 0.57360818672118 | -0.05564580586828204 | 0.0008344556915667738 | 0.3197905661091963 | -0.0639275043333333 | -0.0474335129999999 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_48.0h_STAGEADJ | asn_frac | 4733 | -0.05160125098744888 | 0.0003831448430290192 | -0.07184824072254173 | 7.497259617814477e-07 | 0.08365032295137938 | -0.029874052636927342 | -0.012622777569189933 | acute |
| GSE151189_DHA-CTRL_Cam3II_R539T_48.0h_STAGEADJ_noSchiz | asn_frac | 3227 | -0.037690403607512475 | 0.03227417900204921 | -0.0766301738610871 | 1.3135915596942323e-05 | 0.09707794988399125 | -0.055004755876250006 | -0.023642113491594723 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_8.0h | asn_frac | 5130 | 0.06609488711968392 | 2.1582822550721704e-06 | 0.007470473494220614 | 0.5926896857635742 | 0.00029922874940064073 | 0.0815255222857143 | 0.0284214791666666 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_8.0h_noSchiz | asn_frac | 3610 | 0.07320975179495516 | 1.0673094616244922e-05 | 0.01645563297313082 | 0.32294070794801866 | 0.004814826760110898 | 0.0512543827499999 | 0.0045017261666666 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_8.0h_STAGEADJ | asn_frac | 4760 | 0.016604894087927723 | 0.2520453865314531 | -0.00965018156761666 | 0.5056455599783343 | 0.07729057575599313 | 0.022292754767611108 | 0.0008080115348647614 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_8.0h_STAGEADJ_noSchiz | asn_frac | 3240 | 0.040345945094873985 | 0.021642720699736204 | -0.008322004496585435 | 0.6358406568685444 | 0.05255618006661221 | 0.029296950529495467 | -0.004567737587232159 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_16.0h | asn_frac | 5158 | 0.0948515018824044 | 8.737386282451415e-12 | 0.02966437317708617 | 0.033136581538796485 | 0.0036191071632689607 | 0.0923929179166667 | 0.03982126954166655 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_16.0h_noSchiz | asn_frac | 3632 | 0.1066480499044878 | 1.1680289583980233e-10 | 0.032308734921377785 | 0.05153968440486411 | 0.061960976715098925 | 0.01580516380555555 | -0.0053067089166666005 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_16.0h_STAGEADJ | asn_frac | 4772 | 0.06398410037015734 | 9.717085696687612e-06 | 0.01823287436347409 | 0.2079236019877318 | 0.0488807473216354 | 0.04480030834522932 | 0.008920294573320445 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_16.0h_STAGEADJ_noSchiz | asn_frac | 3246 | 0.09104242643643881 | 2.0382249788551295e-07 | 0.018599320893516777 | 0.28943803985236527 | 0.1368742040684927 | -0.002749499364630613 | -0.020157087041437663 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_24.0h | asn_frac | 5129 | 0.09272828304372488 | 2.859210769814888e-11 | -0.0008336934257744584 | 0.9524008554519228 | 3.5211174918803972e-06 | 0.1408615718333333 | 0.0150828558333333 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_24.0h_noSchiz | asn_frac | 3610 | 0.09894901535243658 | 2.557387274989977e-09 | -0.0021615324882775003 | 0.8967032002561297 | 0.0020232110702927647 | 0.0109579123333333 | -0.0878166458333333 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_24.0h_STAGEADJ | asn_frac | 4753 | 0.07775096322035412 | 8.004682296284512e-08 | 0.004618925996881305 | 0.7502149415633174 | 0.0002723401832681518 | 0.08512272891083722 | -0.0015114584230849458 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_24.0h_STAGEADJ_noSchiz | asn_frac | 3234 | 0.11638614297907425 | 3.164803897067743e-11 | 0.009531331030742139 | 0.5879335143847646 | 0.0016725546558852002 | 0.025772017386573233 | -0.08682018463686479 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_32.0h | asn_frac | 4935 | 0.04325001838754911 | 0.0023740675775421336 | -0.005685117328909357 | 0.6896876251595903 | 0.006493094527406179 | 0.07804907924999996 | -0.0020464405 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_32.0h_noSchiz | asn_frac | 3466 | 0.0364650515678386 | 0.031814208611220175 | -0.02038702465444871 | 0.230165808258005 | 0.09809153941362166 | -0.0101440648333333 | -0.0497514318333333 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_32.0h_STAGEADJ | asn_frac | 4604 | 0.05254713881354487 | 0.0003611387336597664 | 0.0016505168328007904 | 0.9108536756104904 | 0.0027966900254202368 | 0.08053132579087276 | -0.0353004077694375 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_32.0h_STAGEADJ_noSchiz | asn_frac | 3135 | 0.059457339043310437 | 0.0008662815147749839 | -0.008643383450761427 | 0.6285502927913751 | 0.061849927453260296 | -0.045974307635360374 | -0.08613846936812986 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_40.0h | asn_frac | 4976 | -0.09131240144137254 | 1.0957543182473274e-10 | -0.0539971484001287 | 0.00013850458058433678 | 2.7646797435494396e-05 | -0.05787972891666665 | 0.06614557599999996 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_40.0h_noSchiz | asn_frac | 3494 | -0.08804825821713523 | 1.8600899742115714e-07 | -0.06666667635372078 | 8.025099272337686e-05 | 0.006238266494300715 | 0.15195598179999995 | 0.2480223158333333 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_40.0h_STAGEADJ | asn_frac | 4647 | -0.02501479695795373 | 0.08818707226190331 | -0.055610520445287914 | 0.00014894393081888154 | 0.31205562719713453 | -0.01778479455880677 | 0.011321116853475949 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_40.0h_STAGEADJ_noSchiz | asn_frac | 3165 | -0.07294269347539503 | 3.9998975526602714e-05 | -0.06748128852845726 | 0.00014515444593603185 | 0.15472435968380258 | 0.10080702265795148 | 0.14323552167113457 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_48.0h | asn_frac | 5149 | -0.11254225034650774 | 5.522041350382787e-16 | -0.06706872440798783 | 1.4582419490981062e-06 | 2.0880197269427452e-07 | -0.0227226901666666 | 0.053359345833333155 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_48.0h_noSchiz | asn_frac | 3626 | -0.09225040375466018 | 2.622577502146174e-08 | -0.048164958588141425 | 0.003719760758435238 | 0.0001437409344207034 | 0.0134098405555555 | 0.0684013818333333 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_48.0h_STAGEADJ | asn_frac | 4768 | -0.12241425090392805 | 2.2041974023151102e-17 | -0.07719833369509808 | 9.440975919041799e-08 | 8.608440826756707e-07 | -0.08386554919690302 | -0.012641629290221647 | acute |
| GSE151189_DHA-CTRL_Cam3II_revWT_48.0h_STAGEADJ_noSchiz | asn_frac | 3245 | -0.11178904396691787 | 1.7070209635851e-10 | -0.06259250892017926 | 0.0003601272903641345 | 0.00015220228618406727 | -0.05327562364747005 | 0.001131914218338137 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_3.0h | asn_frac | 5174 | 0.0025410503273364995 | 0.8550057884399074 | -0.042501536001860254 | 0.002229678566251145 | 1.427140664264026e-05 | 0.02033690529166665 | -0.0201295874999998 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_3.0h_noSchiz | asn_frac | 3648 | 0.006895037588620198 | 0.6771807139867874 | -0.03807490835846283 | 0.02146410160148567 | 0.013440805668888209 | 0.0238664423333333 | -0.0119091739999999 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_3.0h_STAGEADJ | asn_frac | 4782 | 0.010493614854556578 | 0.4681546734685307 | -0.04107670010848063 | 0.004497334961772099 | 3.0540218584755446e-07 | 0.02983704703067318 | -0.014285898851432909 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_3.0h_STAGEADJ_noSchiz | asn_frac | 3256 | 0.007191670648419298 | 0.6816487392231044 | -0.037280814412966996 | 0.033402058290273345 | 0.004839822425838737 | 0.027729883801261648 | -0.009427168071442745 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_6.0h | asn_frac | 5173 | 0.09300100973488763 | 2.0575977119353147e-11 | 0.05217929755022092 | 0.000173637836878497 | 0.031491106177636896 | 0.03024003511111105 | 0.003569763 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_6.0h_noSchiz | asn_frac | 3647 | 0.09351725085273827 | 1.5313788904852425e-08 | 0.047401934804852636 | 0.004193198176832902 | 0.3616978508226385 | 0.010863461095238 | -0.0041968546666666 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_6.0h_STAGEADJ | asn_frac | 4781 | 0.09051376345874118 | 3.6153512958117405e-10 | 0.05386408983315643 | 0.00019442550990975886 | 0.10558045315589448 | 0.015276877898568119 | 0.0025744382236446107 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_6.0h_STAGEADJ_noSchiz | asn_frac | 3255 | 0.0950476297293053 | 5.543447137070518e-08 | 0.04991520972991192 | 0.0043929353779855245 | 0.5551150906141268 | 0.005031734337455148 | -0.0028836742303127417 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_24.0h | asn_frac | 5172 | 0.08736170207497455 | 3.108384005702964e-10 | 0.03532901571865581 | 0.011056046481180347 | 3.953302618122609e-05 | 0.0995156368333333 | 0.0330790868333333 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_24.0h_noSchiz | asn_frac | 3646 | 0.08318551207574795 | 4.906990770037133e-07 | 0.026769712434750792 | 0.1060624820232517 | 0.013562875648258579 | 0.0288434316 | -0.0539756951666667 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_24.0h_STAGEADJ | asn_frac | 4780 | 0.0950637833519952 | 4.5281996173885e-11 | 0.04355013091378144 | 0.0025988301626002994 | 1.5557636885098815e-05 | 0.09064396961333626 | 0.008193720733233496 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_24.0h_STAGEADJ_noSchiz | asn_frac | 3254 | 0.10724047694953462 | 8.637230506835636e-10 | 0.0452014256260822 | 0.009914692760954988 | 0.003457144615124944 | 0.0037804895859283247 | -0.08764817124944271 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_32.0h | asn_frac | 5167 | 0.06821218655004861 | 9.212795748105397e-07 | 0.05789977668005191 | 3.121318946719724e-05 | 0.020111009740401607 | 0.1098302881111111 | 0.01031035274999995 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_32.0h_noSchiz | asn_frac | 3641 | 0.07408978665211628 | 7.632339595576504e-06 | 0.04658854996881556 | 0.0049273221367419455 | 0.08864525981657433 | 0.1394288033555556 | 0.0434156449166666 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_32.0h_STAGEADJ | asn_frac | 4777 | 0.12590253825999637 | 2.4504373895881083e-18 | 0.07430634485983881 | 2.7253845317385386e-07 | 2.68956427874461e-06 | 0.10783352569607182 | -0.03903307365446812 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_32.0h_STAGEADJ_noSchiz | asn_frac | 3251 | 0.10866956271375885 | 5.228370625076548e-10 | 0.06760026609071218 | 0.00011463165336302434 | 0.010568851224057174 | 0.0639775780747773 | -0.05647016605874012 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_48.0h | asn_frac | 5173 | -0.017390534930207747 | 0.21108714674502804 | 0.021452059978833142 | 0.12290033093271611 | 0.01129134798605792 | -0.0282243877777777 | 0.0252804765 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_48.0h_noSchiz | asn_frac | 3647 | -0.004947240155244065 | 0.7651955081857631 | 0.05193189980268455 | 0.001705523563754694 | 0.036235905634560944 | -0.083117225 | -0.03435754516666665 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_48.0h_STAGEADJ | asn_frac | 4782 | -0.09387328955944878 | 7.81102235437229e-11 | 0.002249782596651995 | 0.8763983770496284 | 5.941107197449183e-07 | -0.1270610903526257 | -0.02798081050867056 | acute |
| GSE151189_DHA-CTRL_Dd2_C580Y_48.0h_STAGEADJ_noSchiz | asn_frac | 3256 | -0.05084810335721553 | 0.0037052138042655103 | 0.0280444334209056 | 0.109608784984999 | 0.0004359191965457981 | -0.15589934575012246 | -0.05681528703659752 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_3.0h | asn_frac | 5173 | 0.06595754123580863 | 2.055132765889657e-06 | 0.07416671772434122 | 9.274072512654801e-08 | 0.021623537987841238 | -0.0034760432499999503 | 0.0318002829999999 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_3.0h_noSchiz | asn_frac | 3647 | 0.09275433200186248 | 2.0040509804604775e-08 | 0.09166361775824199 | 2.9331319772316282e-08 | 0.19892939768631834 | 0.0012465184999999 | 0.033029562416666644 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_3.0h_STAGEADJ | asn_frac | 4781 | 0.07487332778832423 | 2.1831951720315575e-07 | 0.07959426143791218 | 3.5714000497475365e-08 | 0.02805275528382777 | -0.015584067596015109 | 0.017517491187738675 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_3.0h_STAGEADJ_noSchiz | asn_frac | 3255 | 0.10772643869193413 | 7.199021436365205e-10 | 0.10212590649664365 | 5.233115195232421e-09 | 0.15667632265507236 | -0.011736532501791686 | 0.019183299964894943 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_6.0h | asn_frac | 5175 | 0.12226796490673389 | 1.0785035279802308e-18 | 0.09296481560064304 | 2.0763291445404323e-11 | 0.019841246165156878 | 0.0309987996666665 | 0.0175484959999999 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_6.0h_noSchiz | asn_frac | 3649 | 0.12984907245509547 | 3.422250011464741e-15 | 0.09691658714402493 | 4.457794545593339e-09 | 0.04041976758201189 | 0.0311775134999999 | 0.017099631583333302 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_6.0h_STAGEADJ | asn_frac | 4783 | 0.12278144049247064 | 1.5748612914997687e-17 | 0.0973684451688554 | 1.4968404491949052e-11 | 0.027880637669057597 | 0.006010649184180369 | -0.0012233185228681213 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_6.0h_STAGEADJ_noSchiz | asn_frac | 3257 | 0.13855131546032504 | 1.9824224335181046e-15 | 0.10410543587596155 | 2.5972589317900335e-09 | 0.055774019552572385 | 0.011676602497918609 | 0.0021202791716008543 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_24.0h | asn_frac | 5174 | 0.12972698305174166 | 7.344006446396778e-21 | 0.0855567800168906 | 7.094731595867615e-10 | 2.5887733940625354e-08 | 0.19518274608333336 | 0.048584895333333294 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_24.0h_noSchiz | asn_frac | 3648 | 0.13187098664845726 | 1.2743632003440134e-15 | 0.08071132239360775 | 1.0546844330748676e-06 | 6.489809058516291e-06 | 0.1419613891111111 | -0.0292976673333333 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_24.0h_STAGEADJ | asn_frac | 4782 | 0.1591962001190878 | 1.6300631034365787e-28 | 0.10668717401913559 | 1.3948227714289125e-13 | 2.6059383682364694e-10 | 0.152017324844299 | -0.016809614680370594 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_24.0h_STAGEADJ_noSchiz | asn_frac | 3256 | 0.1823861630873365 | 9.469589999850127e-26 | 0.11972611762665453 | 7.192449258515966e-12 | 2.726725891186443e-07 | 0.08170490453669968 | -0.11586483012253504 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_32.0h | asn_frac | 5174 | -0.009628975616550834 | 0.4886447654197412 | 0.0760491145249386 | 4.3293181827304905e-08 | 1.5811391838126854e-05 | -0.02518434424999995 | 0.1051832765 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_32.0h_noSchiz | asn_frac | 3648 | 0.028518395828152545 | 0.0850272659168487 | 0.08874941140481178 | 7.912751463263933e-08 | 0.004460595126169433 | 0.1088471259444443 | 0.1938070380000001 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_32.0h_STAGEADJ | asn_frac | 4782 | 0.0909057305564252 | 3.020694849773886e-10 | 0.0958196123758522 | 3.143908243329787e-11 | 0.9977730191247635 | -0.04007250228966014 | -0.04188139380261503 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_32.0h_STAGEADJ_noSchiz | asn_frac | 3256 | 0.06592535587824291 | 0.00016692124770114667 | 0.10442341947056828 | 2.3351673899172237e-09 | 0.1595568065140327 | -0.031501728732217786 | 0.005894933261761062 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_48.0h | asn_frac | 5174 | -0.045683261613160234 | 0.0010128137078360834 | 0.03161077667884508 | 0.022978106167440088 | 1.3307172272482554e-05 | -0.0594764984166666 | 0.0365972348333333 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_48.0h_noSchiz | asn_frac | 3648 | -0.03262895874678405 | 0.048770199066482216 | 0.06588297139345262 | 6.827817634053918e-05 | 7.929969436249113e-05 | -0.1319208962222222 | -0.0245070583333332 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_48.0h_STAGEADJ | asn_frac | 4782 | -0.12111048355194523 | 4.3228411904078206e-17 | 0.017453251421212526 | 0.2275462600171223 | 5.354121807918912e-12 | -0.19039732135813353 | -0.0376175992763484 | acute |
| GSE151189_DHA-CTRL_Dd2_R539T_48.0h_STAGEADJ_noSchiz | asn_frac | 3256 | -0.06868871334574014 | 8.759833893925744e-05 | 0.052363166802952145 | 0.0028005960336682898 | 4.4964021375767936e-07 | -0.2117399179103886 | -0.07628397846591803 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_3.0h | asn_frac | 5176 | 0.16463380034214742 | 8.974611849748178e-33 | 0.04661186216683716 | 0.0007951521945399122 | 1.0791458765292955e-21 | 0.0419753704583333 | -0.0299274805416666 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_3.0h_noSchiz | asn_frac | 3649 | 0.1791578578000141 | 1.066845365793883e-27 | 0.05196600193218217 | 0.001688590215036192 | 1.1348725755087557e-15 | 0.0486562606666667 | -0.031903542166666646 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_3.0h_STAGEADJ | asn_frac | 4784 | 0.16251413337043294 | 1.1331912746622402e-29 | 0.050906064559039704 | 0.00042772898020677875 | 2.6484519409043113e-19 | 0.06034223309735305 | -0.012937291427446148 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_3.0h_STAGEADJ_noSchiz | asn_frac | 3257 | 0.18201342898386189 | 1.172440889419302e-25 | 0.05802419921644179 | 0.0009231122708008142 | 2.424759568279533e-14 | 0.0695376681313119 | -0.012327731781094696 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_6.0h | asn_frac | 5173 | 0.015411059581351746 | 0.26776838272931197 | -0.004203163994857363 | 0.7624730470815245 | 0.526872925235487 | 0.01405414408333335 | 0.023277378 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_6.0h_noSchiz | asn_frac | 3646 | 0.027681286657532094 | 0.09468183214740397 | 0.002899825344550366 | 0.8610498378613497 | 0.8351949486053113 | 0.0071886958333331 | 0.0140201913333332 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_6.0h_STAGEADJ | asn_frac | 4781 | 0.006398082013664583 | 0.6582845806341463 | -0.006476620466082595 | 0.6543610126489676 | 0.5712005543854792 | -0.011338147889451729 | -0.003977882639437106 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_6.0h_STAGEADJ_noSchiz | asn_frac | 3254 | 0.024860686097024645 | 0.1562423541279527 | 0.0019487617433453413 | 0.9115194024293365 | 0.9094278528278357 | -0.031247643744777773 | -0.010522405732603289 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_24.0h | asn_frac | 5176 | 0.09682826463667803 | 2.9312961642635766e-12 | 0.029789645334672873 | 0.032100509112930165 | 4.283262595370701e-08 | 0.23837499905555548 | 0.0457583263333333 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_24.0h_noSchiz | asn_frac | 3649 | 0.09140112312372563 | 3.18497712716114e-08 | 0.019562302545848738 | 0.23744097114040924 | 0.00014429455465388014 | 0.1467667721666666 | -0.04229646849999995 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_24.0h_STAGEADJ | asn_frac | 4784 | 0.1159655272201387 | 8.555263964285308e-16 | 0.04443481839720275 | 0.0021112835850983784 | 2.5167303611142683e-09 | 0.21738977356680006 | 0.004400616913174761 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_24.0h_STAGEADJ_noSchiz | asn_frac | 3257 | 0.12611328021539334 | 5.070537784160758e-13 | 0.04649504992897528 | 0.007956850465745992 | 3.1001166184252154e-05 | 0.08319238181077673 | -0.10536606465665806 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_32.0h | asn_frac | 5175 | 0.03043637611534186 | 0.028560934538711967 | 0.02587599042632586 | 0.06269966137723072 | 0.1892769100189693 | 0.1444870628492063 | 0.0857170026666666 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_32.0h_noSchiz | asn_frac | 3649 | 0.04554765174505382 | 0.00592541648668915 | 0.015363133021836468 | 0.3535233354799968 | 0.1034141680288957 | 0.3594340316666666 | 0.21187685133333328 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_32.0h_STAGEADJ | asn_frac | 4783 | 0.11223059642671827 | 7.006858773564325e-15 | 0.043702971885863376 | 0.0025018424738959474 | 1.4325493117772424e-06 | 0.2002022602705006 | -0.03199324107335069 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_32.0h_STAGEADJ_noSchiz | asn_frac | 3257 | 0.08011047523776559 | 4.705237055689227e-06 | 0.031047743401692772 | 0.07645338035064564 | 0.0037504426467123695 | 0.2355722328273056 | 0.031807180792376155 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_48.0h | asn_frac | 5176 | 0.02298578676224682 | 0.09822544537727745 | 0.026347481815152046 | 0.05803604391103952 | 0.3477154399447331 | 0.0878187207777778 | 0.0700219871666666 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_48.0h_noSchiz | asn_frac | 3649 | 0.04019489603301313 | 0.015173948900535873 | 0.04979975400056917 | 0.0026203339299714985 | 0.6195350941289723 | 0.0451637157777777 | 0.0442086144166666 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_48.0h_STAGEADJ | asn_frac | 4784 | -0.017326937697959902 | 0.23083145464343366 | 0.01041093553094696 | 0.4715748979016081 | 0.6797806010036456 | -0.042889676843753846 | -0.026771310859481638 | acute |
| GSE151189_DHA-CTRL_Dd2_WT_48.0h_STAGEADJ_noSchiz | asn_frac | 3257 | 0.012248538009624521 | 0.4846872979505592 | 0.03127652964224439 | 0.07430883660715791 | 0.5452591372005784 | -0.06390143282915506 | -0.0372473713252082 | acute |
| GSE225340_DORM_LATE_vs_Rings | asn_frac | 3493 | 0.08286550896948616 | 9.382422987503166e-07 | 0.03183373495489284 | 0.059940874173569365 | 0.8111424759790453 | 0.4858588715794255 | 0.506059586573774 | acute |
| GSE225340_DORM_LATE_vs_Rings_noSchiz | asn_frac | 2542 | 0.0971112754864374 | 9.325141410903682e-07 | 0.030511635163695566 | 0.1240621400406099 | 0.16701666010927807 | 0.5331462768589192 | 0.412389915706739 | acute |
| GSE225340_DORM_EARLY_vs_Rings | asn_frac | 3307 | 0.21651255530071273 | 2.2250025002230146e-36 | 0.11182823357389034 | 1.1284558292093749e-10 | 0.0017792962553496898 | 0.09333780187961116 | -0.06873536094937716 | acute |
| GSE225340_DORM_EARLY_vs_Rings_noSchiz | asn_frac | 2386 | 0.24376367792290066 | 1.2904978013887969e-33 | 0.12132156226204675 | 2.764164453555147e-09 | 0.0006080677207478099 | 0.09185390001242588 | -0.17311597572998894 | acute |
| scRNA_WT_2h_DHAvsDMSO | asn_frac | 5326 | 0.19210477836239728 | 1.9111042661272068e-45 | 0.11293234420891071 | 1.38225457418725e-16 | 2.0239207273867783e-18 | 0.27574074751431876 | 0.17292170024810538 | acute |
| scRNA_WT_2h_DHAvsDMSO_noSchiz | asn_frac | 3768 | 0.17342650963504158 | 7.900944052058101e-27 | 0.10800879153539444 | 2.985804951353479e-11 | 9.049480500347286e-08 | 0.22857604975924062 | 0.1761743488028813 | acute |
| scRNA_580Y_2h_DHAvsDMSO | asn_frac | 5326 | 0.1633358373900055 | 3.639107918984337e-33 | 0.07239505757366693 | 1.2295188878630936e-07 | 1.5810432893168302e-18 | 0.32905462892325854 | 0.17541373754560574 | acute |
| scRNA_580Y_2h_DHAvsDMSO_noSchiz | asn_frac | 3768 | 0.15278463951878918 | 4.0680487945620344e-21 | 0.06130303922880191 | 0.00016633523097972562 | 2.059377898143312e-12 | 0.31628016309357676 | 0.17316283992154968 | acute |
| scRNA_DiD_580YxDHA_2h | asn_frac | 5326 | 0.06652064062755614 | 1.1805023099334443e-06 | 0.019913471608120663 | 0.14620233623058676 | 0.001030844707271204 | 0.059388422697678855 | 0.00844741190246534 | acute |
| scRNA_DiD_580YxDHA_2h_noSchiz | asn_frac | 3768 | 0.07610140881670996 | 2.916071519768864e-06 | 0.015912273781198492 | 0.32881806521222456 | 8.273118606568154e-05 | 0.07504010240343373 | -2.6696278704108067e-05 | acute |
| scRNA_WT_4h_DHAvsDMSO | asn_frac | 5326 | 0.21469521346157353 | 1.4096784369856233e-56 | 0.19182055972314804 | 2.5876097429346313e-45 | 7.808951603610583e-40 | 0.3012469679884937 | 0.03589299376048949 | acute |
| scRNA_WT_4h_DHAvsDMSO_noSchiz | asn_frac | 3768 | 0.19414791164642453 | 2.524571598393272e-33 | 0.18541574873707156 | 1.7126374415462552e-30 | 3.635047497253428e-22 | 0.2140676682610172 | 0.013512949437711974 | acute |
| scRNA_580Y_4h_DHAvsDMSO | asn_frac | 5326 | 0.11756653213889917 | 7.448117487613805e-18 | 0.022627471629493755 | 0.09870505251339107 | 1.2680517365489561e-11 | 0.32955205746098226 | 0.16872115385519137 | acute |
| scRNA_580Y_4h_DHAvsDMSO_noSchiz | asn_frac | 3768 | 0.09484448494531517 | 5.438782883049526e-09 | 0.00852454198817519 | 0.6009000118181809 | 1.1218565661144895e-07 | 0.2401442301225778 | 0.07329916521048663 | acute |
| scRNA_DiD_580YxDHA_4h | asn_frac | 5326 | 0.021065209346454006 | 0.12425974520921072 | -0.06291256752081632 | 4.331672108743528e-06 | 0.5155920226820834 | 0.04201735281899843 | 0.076187770025677 | acute |
| scRNA_DiD_580YxDHA_4h_noSchiz | asn_frac | 3768 | 0.02166221384596605 | 0.18370649777757378 | -0.06364809996986988 | 9.243171721469865e-05 | 0.9558441246574482 | 0.023199672824294026 | 0.018466634525206693 | acute |
| scRNA_WT_6h_DHAvsDMSO | asn_frac | 5326 | 0.0811619117540086 | 3.000909691571412e-09 | 0.0407638616013621 | 0.0029254704200047184 | 6.97471011098418e-16 | 0.5833229840437113 | 0.39781716569285397 | acute |
| scRNA_WT_6h_DHAvsDMSO_noSchiz | asn_frac | 3768 | 0.06322740974984724 | 0.00010285679129912675 | 0.033570785725939405 | 0.039340509650817354 | 1.7415141865758682e-07 | 0.49290533844366635 | 0.3660390693341733 | acute |
| scRNA_580Y_6h_DHAvsDMSO | asn_frac | 5326 | 0.06566665285338145 | 1.615679358589771e-06 | 0.015645991400948803 | 0.25360464980759534 | 3.528719031331708e-08 | 0.3420239283545188 | 0.1067670681029016 | acute |
| scRNA_580Y_6h_DHAvsDMSO_noSchiz | asn_frac | 3768 | 0.04200502194814213 | 0.009916701157600922 | -0.0017156529370807108 | 0.9161545539531479 | 8.956549156758544e-05 | 0.18571397147925417 | 0.005339993491971562 | acute |
| scRNA_DiD_580YxDHA_6h | asn_frac | 5326 | 0.03707222101513254 | 0.00681389845770322 | 0.010018138055367463 | 0.46480057553865795 | 0.1963326001022242 | -0.28501567685728624 | -0.2901555389786541 | acute |
| scRNA_DiD_580YxDHA_6h_noSchiz | asn_frac | 3768 | 0.023814881433979168 | 0.1438574888053075 | -0.003467767690058173 | 0.8314864013176793 | 0.1965856368876724 | -0.31920955370237714 | -0.3611119198294528 | acute |
| Mok2021_Rings_E1_C580Y | asn_frac | 2549 | -0.1034493084729336 | 1.6552477862704026e-07 | -0.0003336854991122774 | 0.9865652766146635 | 0.01849820998096322 | -0.0394999999999994 | 0.06899999999999945 | chronic |
| Mok2021_Rings_E1_R539T | asn_frac | 2605 | -0.1155802596060656 | 3.2960668737321786e-09 | -0.0017044661494960874 | 0.9307089508732822 | 0.0019278421610062497 | -0.1169999999999991 | -0.062499999999999944 | chronic |
| Mok2021_Rings_E2_C580Y | asn_frac | 2191 | 0.031097469966374934 | 0.14563180404706505 | 0.04147504839403063 | 0.05224659340760876 | 0.038607868655223684 | 0.26575000000000015 | 0.1734999999999997 | chronic |
| Mok2021_Rings_E2_R539T | asn_frac | 2205 | -0.057435752473423876 | 0.0069813420873926935 | 0.008987235280458397 | 0.6731807395388911 | 0.279509100308864 | 0.0067000000000003 | 0.0289999999999999 | chronic |
| Mok2021_Troph_E1_C580Y | asn_frac | 2381 | -0.13067377013721096 | 1.5496677320167736e-10 | -0.10357966607677851 | 4.0798098650927523e-07 | 7.493489654790453e-08 | -0.0999999999999996 | -0.0449999999999999 | chronic |
| Mok2021_Troph_E1_R539T | asn_frac | 2387 | -0.11870463011125672 | 5.988308807255937e-09 | -0.0663756927609887 | 0.0011753608522719024 | 1.5756006468295992e-07 | -0.1495000000000006 | -0.0769999999999999 | chronic |
| Mok2021_Troph_E2_R539T | asn_frac | 2034 | -0.05711711480784823 | 0.00998075099368009 | 0.008890500486307134 | 0.6886247937828731 | 0.7497915631284315 | 0.01599999999999905 | 0.030500000000000853 | chronic |
| GSE59099_RvsS_preACT | asn_frac | 5181 | 0.0065733976787392765 | 0.6361856333641907 | 0.016458810084910062 | 0.23622083321836929 | 0.2457539587306995 | -0.015865661499999906 | -0.013204758874999964 | chronic |
| GSE59099_RvsS_preACT_noSchiz | asn_frac | 3652 | 0.0031617287540501478 | 0.8485224875637961 | 0.011433292484250963 | 0.48974135605741276 | 0.6586700589387361 | -0.0059246057500000615 | -0.006747986375000015 | chronic |
| GSE59099_RvsS_Mekong | asn_frac | 5181 | 0.0482711960399796 | 0.0005095066256583635 | 0.049224272159530805 | 0.00039348443124056603 | 0.792378108886556 | -0.015065660500000133 | -0.010812593499999967 | chronic |
| GSE59099_RvsS_Mekong_noSchiz | asn_frac | 3652 | 0.04939544023669971 | 0.00282786578107294 | 0.04409989664080388 | 0.0076892314270285665 | 0.5849844074141443 | 0.0018864883750000297 | -0.0023924224999999633 | chronic |
| GSE226632_dTEinteraction | asn_frac | 5282 | 0.23028800115591236 | 1.589989971052683e-64 | 0.12639330157055595 | 2.950067196811097e-20 | 1.020380764277569e-15 | 0.09801826919241412 | -0.058521296348396194 | other |
