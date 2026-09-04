# 文献前沿与创新边界图

## 1. 本文档用途

Agent 不应把“检索到相关论文”视为完成文献工作。每篇锚点文献必须回答三件事：

1. 它已经证明了什么；
2. 它没有证明什么；
3. 本项目如何在不重复的情况下检验更高一级问题。

详细参考文献表位于 `manifests/literature_anchor.tsv`；原方案书含完整 50 篇参考文献。

## 2. 高水平期刊当前关切

| 期刊/层级 | 疟原虫研究常见主问题 | 证据风格 | 对本项目的约束 |
|---|---|---|---|
| Nature | 不可替代的历史/体内/传播资源，领域级干预 | 独特数据 + 广泛结论 + 体内/公共卫生价值 | 单纯比较基因组不够；只有跨领域原则或体内干预才讨论 CNS |
| Science | 全基因组功能、系统药物进化、生命周期图谱 | 大规模资源服务机制与预测 | “规模”必须形成新功能能力，不能只扩大物种数 |
| Cell | 宿主界面、抗原/抗体与分子机制 | 结构、功能与人体/宿主意义 | 若项目转向侵袭/免疫，需强分子功能闭环 |
| Nature Microbiology | 寄生虫特异调控与可推广微生物学原则 | 跨层组学 + 遗传/生化/自然扰动 | 本项目目标风格：组成历史如何预设翻译应激空间 |
| Nature Communications | 独特资源 + 聚焦功能桥接 | 严密比较/群体数据 + 一项功能验证 | 是机制不够深但资源完整时的现实高位备选 |

## 3. 直接科学锚点

### A. 组成与突变

- **Gardner et al., Nature 2002, 10.1038/nature01097**  
  建立 P. falciparum 极端 AT-rich 基因组事实。未解释跨物种组成转换及功能后果。

- **Hamilton et al., NAR 2017, 10.1093/nar/gkw1259**  
  mutation accumulation 直接支持 P. falciparum G:C→A:T 偏倚，并显示 AT-rich 重复的高 indel。未证明其他物种拥有同一突变机制，也未闭合翻译/适应后果。

- **Otto et al., Nature Communications 2014, 10.1038/ncomms5754；Sundararaman et al. 2016, 10.1038/ncomms11078**  
  说明 Laverania 比较可揭示宿主转换和基因组事件。未系统量化碱基偏倚向蛋白层的机械传导。

**本项目增量：** 先重建多状态组成演化，再以约束反事实分离“背景必然”与“功能残差”。

### B. 密码子、tRNA 与营养

- **Ng et al., Molecular Systems Biology 2018, 10.15252/msb.20178009**  
  将 tRNA epitranscriptomics、codon bias 与蛋白表达相连。主要是单物种机制。

- **Li et al., iScience 2024, 10.1016/j.isci.2024.111167**  
  将氨基酸需求、tRNA 供应、血红蛋白营养与增殖基因表达相连，并做有限跨物种 proof-of-concept。未控制完整系统发育、同源位点和组成反事实。

- **Small-Saunders et al., Nature Microbiology 2024, 10.1038/s41564-024-01664-3**  
  证明 tRNA 修饰重编程和 Lys codon-biased translation 参与青蒿素存活，并有 PfMnmA 扰动。已经占据“动态 tRNA 修饰导致耐药”这一具体机制。

**本项目不能重复：** “耐药株会调节 tRNA 修饰”或“P. falciparum 有密码子偏好”。  
**本项目拟回答：** 长期谱系组成演化是否预先限制/扩展了不同物种可调用的翻译应激程序。

### C. 生命周期与翻译抑制

- **Howick et al., Science 2019, 10.1126/science.aaw2619；Dogga et al., Science 2024, 10.1126/science.adj4088**  
  提供完整生命周期/性发育单细胞表达框架。

- **Lindner et al., Nature Communications 2019, 10.1038/s41467-019-12936-6**  
  孢子体成熟中存在两波翻译抑制，证明阶段间 RNA 与蛋白输出可以解耦。

- **Afriat et al., Nature 2022, 10.1038/s41586-022-05480-z；Yan et al., Nature 2025, 10.1038/s41586-025-09653-0**  
  说明肝期和蚊期仍是前沿盲区，且高水平工作要求空间/宿主交互或体内验证。

**本项目增量：** 将序列组成残差作为阶段翻译输出差异的先验，而不是重新做阶段聚类。

### D. 功能基因组与耐药

- **Zhang et al., Science 2018, 10.1126/science.aap7847**  
  P. falciparum 饱和突变必需性基线。

- **Luth et al., Science 2024, 10.1126/science.adk9893**  
  724 个药物选择基因组和 118 个化合物形成系统耐药图谱。可用于校准突变可达性和候选逃逸，而不是让本项目再做泛化 drug-target list。

- **Elsworth et al., Science 2025, 10.1126/science.adq6241；Oberstaller et al. 2025, 10.1126/science.adq7347**  
  P. knowlesi 必需基因组及跨物种 essentiality rewiring，直接说明 P. falciparum 机制不可默认普遍。

- **Birnbaum et al., Science 2020, 10.1126/science.aax4735**  
  K13 定义的内吞路径连接青蒿素耐药，显示耐药是细胞状态、代谢与蛋白稳态问题。

**本项目增量：** 解释哪些长期编码架构与现有必需性/药物进化结果一致，并提出跨谱系适用边界。

### E. 群体与转化

- **Ibrahim et al., Nature Communications 2024, 10.1038/s41467-024-55102-3**  
  P. malariae 251 个基因组与 P. knowlesi ortholog replacement 验证药敏，是“跨物种发现—可培养物种验证”的重要先例。

- **MalariaGEN Pf8, DOI 10.12688/wellcomeopenres.24031.1；Pv4, DOI 10.12688/wellcomeopenres.17795.1**  
  群体数据足以支持自然变异、选择与 CNV 验证。不能以“寄生虫没有大数据”为前提。

- **Billows et al., Nature Communications 2026, 10.1038/s41467-026-73006-2**  
  全球尺度群体适应进一步抬高了“仅描述群体结构”的门槛。

## 4. 学界最关心的接口

本项目只有连接到以下至少一个接口，才有 Nature Microbiology 级意义：
- 青蒿素耐受、恢复/持留和蛋白稳态；
- 不同物种/谱系的功能可迁移性；
- 性发育、孢子体、肝期或传播中的翻译切换；
- 基因组可塑性导致药物、诊断或免疫逃逸；
- 寄生虫特异 tRNA/翻译机制与宿主选择性。

## 5. 创新性审计问题

每个主要结果在进入主图前必须回答：
1. 是否只是“高 AT → AT-rich codon/AA”的重述？
2. 是否已被 iScience Figure 6 或 Nature Microbiology 2024 直接完成？
3. 是否有系统发育独立性？
4. 是否控制了 LCR 和注释？
5. 是否预测了一个独立功能表型？
6. 是否产生一个既往工作无法给出的边界或可证伪预测？

任一结果仅满足前两层描述时，进入 Extended Data 或资源表，不占据主标题。

## 6. 文献更新制度

- M0：完成系统检索与现有工作矩阵；
- Gate B：对选定主轴进行一次针对性更新；
- M4：在锁定稿件故事前再次检索最近 12 个月；
- 发现高度重合研究时立即创建 change request，不等待写稿阶段。
