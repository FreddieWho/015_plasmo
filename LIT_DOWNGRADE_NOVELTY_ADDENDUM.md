# 补遗：降级结论 5 篇文献复核（D-015，2026-09-04）

对应 `LIT_DOWNGRADE_NOVELTY_CHECK.md`（D-014，当时无检索工具，5 篇标存疑）。
本次用更新后搜索工具全文核验（Europe PMC/PMC 全文 + M3-02 本地 PDF 抽取），存疑全部解除。

## 逐篇边界（已验全文）

- **Hamilton 2017 NAR (PMC5389722)**：单物种。6 个 Pf 分离株、4 年、279 克隆的 MA，无跨物种分解。indel 富集 AT 重复/非编码区，外显子内保框（3 的倍数）富集。我方引用边界正确；跨物种 LCR 因果检验无人做过。
- **Otto 2014 Nat Commun (PMC4166903)**：HKA/MK/Ka/Ks 基因级宿主适应扫描。无密码子偏好分析，无组成方差分解。
- **Sundararaman 2016 Nat Commun (PMC4804174)**：密码子仅用于建树比对（三位点去除）；LCR 屏蔽（tblastx ≥20aa/95% + segmasker）仅为建树比对 QC，非保留/屏蔽因果双轨。我方双轨仍新。
- **Li 2024 iScience (PMC11544085)**：Pf 内 TPM 加权 demand + tRNA supply S/D + wobble + anota2seq + 放线菌素 D decay——与 M3-01 的 Pf 内方法同构，**必须引用、防重复宣称**。跨物种仅 Fig6（糖酵解/PPP 的 AA 组成 vs 宿主 HB 相关，Mann-Whitney，无系统发育控制）+ Fig3E（top5% Ile/Leu-rich 基因的 ortholog 组 SD 检验，43/53 Ile-rich >3SD）。全蛋白组 CLR + LCR 双轨 + LOO + GC3/dinucleotide null 比较仍是我方独有。
- **Small-Saunders 2024 Nat Microbiol (PMC11153160)**：单物种 Dd2（M3-02 已从本地 PDF 抽取验证），无需复核。

## 事后敏感性分析（post-hoc，不改 Gate 判据）

Pf-vs-Pk 残差向量相关性（M3-01 主表，36 家族）：
codon 层 Pearson 0.879 / Spearman 0.430（Pearson 由 L_CTN/L_TTR/S_AGY/R_CGN 四个拆开亚家族的大点驱动）；
AA 层 Pearson −0.588 / Spearman −0.516（镜像=组成轴形状，AT 富 vs GC 富反向拉，强化降级结论）。
|dresid|>0.2 的一致项仅 4 个拆开亚家族；K 残差 +0.0326/−0.0073，量级≈零。
解读：“共享程序”检验即使阳性也与谱系特异 H3 矛盾，救不了 H3 跨物种版；AA 负相关反而支持组成轴解释。

## Figure 2 最终定位（不变）

正文校正位：效应量越大越好，只讲“背景扣干净了”，不讲适应性发现；(3) LCR/模型比较放 Extended。
