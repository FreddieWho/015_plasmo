# M5 Convergence Report — M4R Round 1 (2026-09-05, 主Agent统一审计)

**问题：** 当前最强故事是 D / D+A / D+C / D+A+C / 仍不足回降级？（不为期刊选答案）
**输入：** M4R-D（worker, CONDITIONAL）/ M4R-A（worker, CONDITIONAL）/ M4R-C（worker×2超时后主Agent收尾, CONDITIONAL）。
**方法纪律：** association≠mechanism；同实验多模态=orthogonal非independent；探索标签保留。

## 1. D：能否成为 main evolutionary backbone？—— 能，CONDITIONAL

- D1 branch-aware（24,520 mapped subs）：informative branches总体concordance 0.583（p=3.2e-83）；coupled 0.791 vs uncoupled 0.628，Fisher OR=2.25 p=9.5e-10；multi-nt>1-nt accessibility梯度一致。**Branch层已强于terminal层**（terminal保留descriptive）。
- D2 architecture诚实mixed但利好general principle：core基因~99%耦合turnover落在**非LCR结构化序列**（convention caller在蛋白上0%双组；homopolymer OR=17.5但counts tiny不承重）+ radical chemistry（Grantham 77.3 vs 66.9）→ churn不是repeat故事，general-principle读法成立。
- D3 matched controls：14/14耐药位点uncoupled；conservation梯度真实（cons=1.0耦合0% vs 0.5处9.1%）→ **保留弱结论**（constrained resistance sites vs turnover elsewhere），不称双轴。
- D5：AP2 MH OR=27.4（p=5.4e-11，2/2层同向）+ PUF OR=3.03 + CHROM null + essentiality约束（AP2 blood-stage可耐受插入→lifecycle-conditional框架）。Counterfactual成立。
- 短板：D6仅pilot（AP2/PUF不在core-181内，结构性局限；4 OG方向一致但underpowered）；D4仅boundary（Pf8未算）；branch-LOO未跑。
- **结论：D可扛Figure 2–3（evolutionary backbone），pending round-2形式化。**

## 2. A：是否超越"regulators are Asn-rich"？—— 未达，CONDITIONAL

- 富集真实：AP2 OR=34.3（13.8–84.9）q=1.5e-16 + transcription_reg/sexual/CCR4-NOT/RNA-binding同向 + 调整OR=2.59；但反事实仅high-LCR层成立、zero-LCR n.s.；严格PUF集不可判；chromatin/proteostasis null。
- Lifecycle-state consequence：GSE75795方向存在但小且BH未过（n=1/性别，探索标签）；ChIP单研究弱阳性（OR=1.55）+ AP2-I反向；GSE222586/220039未入库；**本轮不收候选shortlist**。
- GCN5 2026阳性对照在手（生长崩溃~15× + Py互补），cite-and-complement合规。
- **结论：A是functional pillar条件性候选；lifecycle co-option仍为hypothesis，需round-2（lifecycle两集+同背景perturbation表达+正式GO/InterPro复现）。**

## 3. C：stage/confounder修正后剩多少？—— 剩一个partition figure，CONDITIONAL

- Stage修正后acute存活：microarray早期/中期（Dd2_WT 3h 0.165→0.163；R539T 6h→0.123；stage_R²~1e-4，stage非驱动）+ scRNA独立同向（2–4h 0.12–0.22，lenadj~0.11）+ persistence LATE第三向（0.083）。效应小（0.06–0.23），length衰减大但不归零。
- Acute vs chronic：transcript两端正 vs Mok蛋白负（−0.10~−0.13）vs GSE59099转录null（0.007 n.s.）→ 分区hypothesis成立，节约在蛋白层。
- DiD基因型交互微弱（0.02–0.07）→ 无K13特异MoA证据；M4R-X01 gap维持；dTE规则2（方向一致agree 0.924，q仅1基因→弱支持保留）。
- 退出边界未触发 → C保留为partition/boundary figure候选，非MoA figure。
- **结论：C提供第二functional context（Extended/Fig5候选），不独立承重。**

## 4. Integrated verdict（八维度简评）

| 维度 | 评估 |
|---|---|
| effect size | D1 OR 2.25 / AP2 OR 34（CI下限13.8）强；A-state/C-acute小（0.06–0.23） |
| phylogenetic robustness | D1 branch层天然phylo-aware；branch-LOO待补；A/C单物种（Pf） |
| confounder robustness | length/LCR/schizont/stage多轨已跑；length衰减是最大诚实减分项 |
| cross-dataset independence | C-acute双transcript+第三方persistence；A-lifecycle缺第二独立集；D缺Pf8 |
| biological specificity | CHROM/proteostasis nulls + AP2-I反向 + chronic转录null = 非万能富集 |
| mechanism discriminability | 未达（四机制/acute-chronic均只到分区描述；translation层缺失） |
| novelty | NT-1/PA-1已fence；branch-aware+partition+axis-separation组合仍新 |
| Nat Micro relevance | D general principle + GCN5 lifecycle proof + ART partition对齐期刊接口 |

**六层充分性：** L1✓ / L2✓+branch / L3✓弱版 / L4✓（structured churn+regulator富集）/ L5半（GCN5单基因座实证+探索性state方向+ART分区）/ L6待（机制判别/最小实验）。

## 5. 五选一

**当前答案：D + A方向（D+A leading），D+C为第二context，D+A+C仅当round-2两表型线同指一机制时采用；降级暂不执行。**

- D backbone CONDITIONAL-strong（Figure 2–3可立）。
- A pillar CONDITIONAL（Figure 4候选，待round-2 state-consequence转正）。
- C partition CONDITIONAL（Figure 5/Extended候选，MOA不指望）。
- T01维持SUPERSEDED_FOR_NOW；Nat Micro重回主轨道但**尚未充分**（L5/L6待补）——诚实状态，非选答案。

## 6. Round-2（收敛前最后bounded增量，不扩主题）

1. D6形式化（≥2独立transitions预注册pole对比；AP2/PUF wide OGs）+ Pf8罕见变异边界 + branch-LOO。
2. A：GSE222586/GSE220039入库映射 + PfAP2-P/PfPuf1 processed表达 + InterPro正式富集复现；成立则收1–3候选。
3. C：维持现状（不再加模态）；仅若acute translation×K13数据出现才升级。
4. Round-2后开M5终审：A-state转正则D+A成稿；A-state失败则D+C partition或诚实降级（T01重激活）。

## 7. 对PLAN假设的影响

- H-D（partition）：round-1支持（structured churn + 弱版axis-separation + AP2/PUF counterfactual）。
- H-A（lifecycle co-option）：未证实也未证伪，富集侧强、state侧弱。
- H-C（ART context）：分区描述成立，机制未立；supply措辞维持封顶。
