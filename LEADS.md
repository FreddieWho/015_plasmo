# LEADS — 扫过未深究、有潜力的方向（只追加，不删）

2026/9/3
## L-001 AT-rich 重复→indel/CNV 可塑性 fallback（H5）
M2 未测重复/断点耦合，若 M3 翻译失败可转主轴。
最小验证：重复背景 vs MalariaGEN CNV 断点富集，控制长度/AT/区室。
状态：待挖掘

2026/9/3
## L-002 孢子体翻译抑制触发 D022/D023（H3 阶段扩展）
M2 残差若富集传播模块才下载 GSE113582 + PXD009726–29。
最小验证：残差模块与孢子体两波抑制基因集重叠。
状态：待挖掘

2026/9/3
## L-003 P01 vs Salvador-I 亚端粒差异（H1 稳健性）
株差可能被误判为物种差异，需双株并跑组成与残差。
最小验证：双株 GC/基因数/残差并排审计图。
状态：待挖掘

2026/9/3
## L-004 全 MSA + 祖先重建升级 M2-02 site 层（H2 精度）
当前 site 层空数值占位，中心星形 + 100 OG 抽样。
最小验证：MAFFT 全量比对 + 简约祖先，限 181 core。
状态：待挖掘

2026/9/3
## L-005 Pf8/Pv4 群体变异验证候选容忍度（H4 逃逸）
自然变异可做低逃逸论证，不做全量 FASTQ。
最小验证：候选位点 VCF/Zarr 等位频率 + CNV 查询。
状态：已并入主线（R3，Pv4 作辅线旁证，D-008）

2026/9/4
## L-006 降级结论 5 篇近年疟原虫文献复核（D-014 后续）
researcher 无检索工具，Hamilton/Otto/Sundararaman/Li/Small-Saunders 是否做过跨物种组成分解或 LCR 因果检验未知。
最小验证：有检索后补 5 篇 DOI+摘要筛查 + LCR/dinucleotide 二次检索，复判 Figure 2 去留。
状态：已并入主线（D-015，全文验完，存疑解除）

2026/9/4
## L-007 AAT/Asn 轴的蛋白层确认（机制链第 3 环）
L1 双阳 + L1c 供应（Asn-tRNA 饥饿下降 6/44）+ L1d 功能（Asn 富集=AP2/PUF 调控层）已成，缺蛋白量层的独立确认。
若 TMT/蛋白组里 Asn 富集蛋白在应激下同向变化，机制链闭环，可重新评估期刊上限。
最小验证：D008 Source xlsx（用户浏览器下载）或 Mok 2021 processed 蛋白组；按 AAT top-decile 分组比较蛋白响应。
状态：待挖掘

2026/9/4
## L-008 PXD047875（Bozdech 2025）热稳定性数据
researcher 列为候选第三数据；实查为 ITDR 热位移设计（8-23GB/包 raw），非丰度时间序列。
对我们的丰度检验是错误设计；但若未来问"Asn 富集蛋白是否更易热失稳"（应激脆弱性另一角度）可用。
最小验证：暂不投入；仅在需要热稳定性轴时重新评估。
状态：已放弃（设计不匹配）

2026/9/5
## L-006 D2b 全蛋白组 architecture 归属（不限 core-181）
core-181 采样框排除 LCR 重/快进化蛋白，"churn 不在 LCR"部分由框架构造。
若 D2b 显示框外 churn 实为 LCR 富集，Figure 3 与 general-principle 强度都要降。
最小验证：全 Pf 蛋白组 composition-coupled 候选 × LCR/IDR/domain 归属表，一次工作单元。
状态：待挖掘（round-2 硬项）

2026/9/5
## L-007 严格 PUF 集构建（PlasmoDB annotation/orthology，非 regex）
round-1 严格 PUF 仅 2 基因不可判；D 线 56 基因宽集口径不能当 PUF 证据。
若严格集（~10–20 基因）富集成立，CLM09 才能写 AP2/PUF。
最小验证：InterPro/ortholog 组 PUF 名单 + 同样富集流程复跑。
状态：待挖掘（round-2 硬项，并入 A）

2026/9/5
## L-008 branch-aware LOO 系统发育敏感性
D1 concordance（OR=2.25）未做 leave-one-clade-out；单枝驱动风险未排除。
最小验证：7 次 LOO 重跑 concordance/Fisher，几小时计算。
状态：待挖掘（投稿前必做）

2026/9/5
## L-009 dTE interaction 正式 anota2seq 复跑
round-1 interaction 为 OLS+t（df 小，q 仅 1 基因）；规则 2 保留弱支持。
anota2seq（作者原方法）若同向存活可升级为正式 interaction evidence。
最小验证：R anota2seq 一次安装 + 12 样本矩阵，半天。
状态：待挖掘（可选，非 gate）
