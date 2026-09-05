# 分析合同与模型边界

## 1. 总原则

方法选择服从科学分解，而不是“现代工具清单”。每个主要分析必须包含：

1. 明确的统计单位；
2. 明确的零模型/反事实；
3. 效应量与不确定性；
4. 系统发育与数据质量控制；
5. 至少一项独立敏感性分析；
6. 结果为阴性时的解释边界。

## 2. 数据质量合同

### 基因组与注释
- 冻结 assembly accession、annotation release、文件 checksum；
- 核、线粒体、顶质体和未定位 contigs 分开；
- 记录 N50、总长、N 含量、基因数、完整 CDS 比例、BUSCO/核心 ortholog 回收率；
- 物种比较优先使用高可信单拷贝 orthogroups；
- 不因一个物种注释更多而解释为真实基因家族扩张；
- 对关键候选至少交叉检查 NCBI 与 PlasmoDB/文献注释。

### ID 与序列
- 所有 gene/transcript/protein ID 映射具有版本与多对一规则；
- 主分析每个基因只使用预先定义的代表转录本；
- 终止密码子、内部 stop、frameshift、过短 CDS 和假基因需显式标记；
- 遗传密码表按核/细胞器分别确认；
- 原始和清洗后序列均保留。

## 3. Aim 1：组成演化

### 必需描述层
- 全基因组、CDS、内含子/UTR/间区（可得时）；
- GC1、GC2、GC3、四倍简并位点；
- 单/二核苷酸上下文、homopolymer、短串联重复；
- 基因组区室和距离端粒/亚端粒；
- LCR 保留与屏蔽。

### 推断边界
- 参考基因组：组成；
- 多物种树：祖先组成与替换；
- mutation accumulation/de novo/低频变异：突变谱；
- 不得通过 neutrality plot 单独宣称突变与选择贡献比例。

### 系统发育
- 主要物种树由高可信单拷贝核心 orthologs 建立或经权威树验证；
- 对关键转换检查基因树冲突和长枝；
- 生态关联必须使用系统发育比较或等价控制；
- 至少执行 leave-one-clade-out。

## 4. Aim 2：同义密码子反事实

### 目标量
估计在给定氨基酸和关键背景后，真实同义密码子使用相对期望的偏离。

### 主零模型至少保持
- 氨基酸序列；
- 基因长度和密码子位置；
- 根据问题选择保持全基因/局部 GC、二核苷酸或相邻密码子背景；
- 物种/谱系结构。

### 必需混杂控制
表达量（仅深挖物种）、基因组区室、低复杂度、转录本选择、tRNA gene/abundance 数据可得性、起始附近序列与重复。

### 传统指标定位
RSCU、ENC、CAI、GC3、PR2、PCA/对应分析用于描述、QC 和旧文献对齐；不能独立支撑适应性。

## 5. Aim 2：氨基酸与蛋白反事实

### 主零模型至少保持
- orthogroup 与可比同源位点；
- 祖先状态或可行替换路径；
- 蛋白结构/功能环境的可用近似；
- LCR 状态与蛋白长度；
- 组成闭合效应。

### 组成型数据
氨基酸频率总和为 1。主要分析应使用 Aitchison/CLR/ILR、multinomial、Dirichlet-multinomial、logistic-normal 或其他明确处理闭合约束的方法；原始比例相关仅作可视化。

### 结构分层
至少区分结构域、无序/LCR、信号肽、跨膜、成熟蛋白；不能把定位肽或膜蛋白的物理化学需求解释为营养适应。

## 6. Aim 3：翻译供需

表达加权需求可抽象为：

\[
D_{c,s,t}=\sum_g E_{g,s,t}N_{g,c}
\]

其中 `E` 是物种 s、阶段 t 的表达权重，`N` 是基因 g 的密码子或氨基酸计数。

**最低要求：**
- 不把转录本丰度直接等同于翻译速率；
- tRNA gene copy、abundance、charging、修饰和 wobble 是不同证据层，不可混为一个“tRNA supply”；
- 不同技术/物种的数据先在物种内标准化，再讨论跨物种；
- 血红蛋白组成只是宿主营养第一近似，不是实测通量；
- 主候选需预测 RNA—蛋白偏离、tRNA/核糖体、mRNA 稳定性或扰动中的至少一项。

## 7. 模型层级

### 一级：主推断
- 描述统计与效应量；
- 系统发育广义最小二乘/混合模型或等价方法；
- 祖先状态和谱系转换；
- 位点/密码子模型；
- 约束随机化和置换；
- 多层线性/广义线性模型。

### 二级：增强
- 贝叶斯层级模型；
- 结构/功能分层；
- 非线性广义加性模型；
- 专门的选择/进化模型。

### 三级：候选排序，不承担主要因果
- 随机森林/梯度提升；
- 蛋白语言模型 embedding 或变异效应评分；
- AlphaFold/口袋/网络排名。

**禁止：** 因一级模型结果不支持假设而直接升级三级模型“寻找信号”。

## 8. 验证单元

- 跨物种泛化：leave-one-species-out 或 leave-one-clade-out；
- 通路/基因集：按 orthogroup/物种层 bootstrap，不把单个密码子当独立重复；
- 候选稳定性：替代注释、LCR 屏蔽、代表转录本规则、物种面板变化；
- 生态趋同：至少两个独立转换，或清晰的近缘匹配对照；
- 多重检验：同时报告 FDR、效应量、置信区间和方向一致性。

## 9. 候选进入机制层的准入

候选至少满足：
1. 残差在主要敏感性分析中方向稳定；
2. 不是由单一物种、蛋白长度或 LCR 驱动；
3. 有第二证据层支持；
4. 与领域核心表型有直接路径；
5. 存在可区分的竞争解释；
6. 不是对既有 Nature Microbiology 结论的简单复述。

## 10. 结果交付合同

每个分析目录必须包含：
- `README.md`：问题、输入、方法、输出、结论边界；
- `params.yaml` 或等价配置；
- `input_manifest.tsv` 与 checksums；
- 可复现脚本/工作流；
- 主结果表，不仅是图片；
- QC 与敏感性结果；
- `claim_impact.md`：支持/削弱哪些 Claim 和 Gate；
- 环境/软件版本。

## 11. M4R 分析增补（2026-09-05，D-036）

- **branch-aware site analysis：** 高质量 ortholog alignment → ancestral AA reconstruction → substitution 映射到 species-tree branches → 每 branch 定义 composition change → 判断 observed AA replacement 是否与 branch GC/AT change 同向；branch 为独立信息单位；不可靠位点标记 uncertainty 后排除，不用复杂 Bayesian machinery 强救。terminal correlation 保留为 descriptive layer。
- **codon accessibility：** 每替换报告最少 nucleotide changes、GC→AT/AT→GC direction、synonymous/nonsynonymous accessibility。
- **structured/IDR/LCR partition：** composition-coupled vs uncoupled 位点映射到 structured domain / conserved domain / IDR / LCR / homopolymer-repeat / linker / TM / signal / catalytic-binding（如可靠）/ exposure-confidence（仅可靠候选）；至少控制 protein length、local conservation、orthogroup、LCR/IDR availability、gene functional class、annotation quality。不得预设信号集中 IDR；structured 内的 radical churn 如实解释。
- **matched adaptive-site controls：** 每个 known adaptive/resistance site 匹配 conservation、protein/domain、structural context、essentiality/constraint（如可得）、alignment quality 相当的 background sites，再比 coupling/turnover/accessibility/Grantham/population variation。若匹配后差异消失，保留弱结论（canonical resistance mutations occur at highly constrained sites, whereas composition-driven turnover dominates a different part of sequence space），不硬称“双适应轴”。
- **whole-protein Asn feature decomposition：** 禁止只用 whole-protein Asn fraction；至少同时测试 whole / structured / IDR / LCR / poly-Asn tract Asn + matched protein controls，判断 ART signal 归属（Asn chemistry vs LCR/IDR vs length vs developmental programme vs general stress group）。
- **stage-adjusted ART：** DHA/ART time-course 能推断 developmental age 时用高分辨率 IDC reference 分离 nominal time 与 inferred stage，重算 stage-adjusted response；不能可靠推断时明确限制。acute 与 chronic 分开分析，禁止把反号直接解释成 mobilization/conservation 机制（仅 hypothesis）。
- **dTE evidence hierarchy：** 现有 GSE226632 ratio-based dTE 保留，状态 DISCOVERY/SUPPORTING；normalized polysome/total、condition×fraction interaction、anota2seq 或等价 replicate-aware model 为 upgrade 而非 hard validity gate。新模型同向存活→升级；同向但 significance 下降→保留弱支持；系统性翻转→调查并降级；无法重建→保留探索标签。不把方法学完美主义变成删除发现的理由。
