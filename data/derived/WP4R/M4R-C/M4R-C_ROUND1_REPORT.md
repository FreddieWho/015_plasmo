# M4R-C ROUND1 REPORT — ART/stress first-gain analysis (2026-09-05, 主Agent收尾)

- task_id: M4R-C-ROUND1；packet: docs/tasks/M4R-C_TASK_PACKET.md
- verdict: **CONDITIONAL**（acute 端多数据集同向弱阳性 + chronic 蛋白反号 + chronic 转录 null；partition figure 条件成立，MOA 不成立）
- outputs: data/derived/WP4R/M4R-C/（README.md below, params.yaml, input_manifest.tsv, M4RC_*.tsv ×9, M4RC_dTE_upgrade.json, claim_impact.md）
- 执行说明：两轮 subagent 连续超时后由主Agent接管完成。根因已定位并修复——scipy 1.11.3 `mannwhitneyu` exact 路径在高度 ties 组（Cam3II 单重复描述性组）组合爆炸；全线强制 `method='asymptotic'` + 退化组跳过（nunique<10）后全管线 ~8 分钟跑通。此前两轮的中间结果（stagevec/stage-adjusted DHA表/MCA伪bulk/protein features）全部复用，未重算。

## 1. 输入版本

- GSE151189 microarray（D010，27 bg×time 组，raw + MCA伪bulk stagevec残差 STAGEADJ + noSchiz 双轨）
- GSE225340 persistence FPKM（M4R-C01：Rings×3 vs D1–4/D5–12）
- K13 scRNA 14 DGE伪bulk（M4R-C03：MRA1250 WT / MRA1251 580Y × Ctrl/2h/4h/6h × DHA/DMSO；transcript-only）
- Mok 2021 MOESM（manual_inbox/M04：7 genotype contrasts；chronic protein 端）
- GSE59099（M4R-C02b：pre-ACT R=378 vs S=361 + Mekong子集；chronic transcript 端；halflife双键合并修复）
- GSE226632 12 htseq counts（condition×fraction interaction；replicate-aware upgrade尝试）
- 方法：spearman(Asn feature, response) + length-adjusted（log10 length残差）+ Asn-top10 MW + noSchiz；BH per contrast（feature_tests_all.tsv）

## 2. 关键效应量（Asn frac；lenadj = length-adjusted；详见 M4RC_feature_decomposition.tsv）

### C1 stage-adjusted DHA（microarray）：部分存活——如实计数 16/27，非全部存活
- 27 组 bg×time 中，raw 显著 21 组；STAGEADJ 后显著且同向 **16/27**（Dd2 背景稳定：WT 3h 0.165→0.163、R539T 6h 0.122→0.123、C580Y 6h 0.093→0.091；Cam3II 混合：16/32h 存活、**24h 与 40h、revWT 8/40/48h、Dd2_R539T 48h 共 6 组原显著转弱或转负**）。
- **弱参考声明：** stagevec 为 MCA 伪bulk 单一 ring-vs-schizont 轴，stage_R² 仅 ~1e-4–3%——"stage 不是驱动者"只在"该弱参考下信号不变"意义上成立，不能排除更强 stage 定义（per-sample inferred age）下的混杂。本轮满足 C1 的"能推断则调整"下限；更强的 stage 解析是 round-2 可选项。
- 40/48h 晚期翻负（恢复/死亡动力学，非急性响应）。
- 衰减项：lenadj后 0.16→0.05（WT 3h）：**length confound 大，但不归零**。
- 归属：asn_out_lcr ≈ asn_in_lcr（0.14 vs 0.17）→ whole-protein Asn + length，不是 LCR-specific。

### C1b scRNA acute（独立数据集，同向）
- WT 2h 0.192 / 4h 0.215 / 6h 0.081；580Y 2h 0.163 / 4h 0.118 / 6h 0.066；lenadj后 2–4h ~0.11存活。
- DiD（genotype×drug）：2h 0.067、4h 0.021 n.s.、6h 0.037 → **药物响应存在，基因型交互微弱**。

### C2 acute vs chronic 分区（禁止反号即机制）
- Acute transcript两端（microarray + scRNA）+ persistence LATE（rho=0.083 p=9e-7，asn_out_lcr驱动）：同向（正）。
- Chronic protein（Mok 5/7负，最强 −0.131 p=1.5e-10，lenadj存活）：反号（负）。
- Chronic transcript（GSE59099 RvsS：0.007 n.s.；Mekong 0.048 p=5e-4 tiny）：**null** → 节约发生在蛋白层，稳态mRNA不可见。重要discriminator。
- 判读：acute动员 vs chronic节约的分区描述成立（hypothesis级），不是已证机制。

### C3 feature分解结论
- Signal归属：whole-protein Asn + protein length；不是 poly-Asn特异（polyN_maxrun与lcr_frac同量级），不是LCR-specific（out≈in）。
- AP2/PUF MW结果见decomposition表is_*行（MW_AsnTop10_p列）；本轮不做regulator富集宣称（归A线）。

### C4 dTE interaction upgrade（4规则执行）
- sign_agreement=0.924，spearman=0.996（p=0）→ 方向与ratio一致；但 q<0.05 n=1 → significance崩。
- 规则2：**同向但significance下降 → 保留弱支持**。ratio-based dTE维持DISCOVERY/SUPPORTING，不升级不删除。

## 3. Verdict理由（CONDITIONAL）
- 强侧：Dd2 背景 stage 修正后存活（弱参考声明见 C1）+ scRNA独立同向（mRNA 层）+ persistence第三向（dormancy，非 acute DHA）+ dTE方向一致 + chronic蛋白反号而转录null（layer-discriminating）。
- 弱侧：效应小（rho 0.06–0.23）+ length衰减大 + STAGEADJ 6/27 组转弱 + DiD微弱 + 无translation层数据（M4R-X01 stands）+ chronic只有蛋白端 + acute 三个数据集全为 mRNA 层（跨数据集独立性评分降半档）。
- 退出边界未触发（stage-corrected非弱/阴；null的独立ART数据集只有GSE59099一个transcript端，不足≥2）→ C继续，但只作partition figure候选。

## 4. 什么能改变verdict
- →STRONG（D+C条件）：acute translation×K13数据（Ribo/polysome/proteome×genotype；现不存在）或第二acute蛋白数据集同向。
- →WEAK/BOUNDARY：length+expression匹配后acute信号消失；或第二独立chronic蛋白数据集null。

## 5. QC/失败记录
- MWU exact-hang根因+修复（见首节；教训：ties重的单重复描述性组禁用exact方法）。
- GSE59099 halflife双键（1043+110互补）按位置合并；R/S定义hl≥5/≤3 pre-ACT。
- dTE interaction：condition×fraction+rep2/3协变量OLS，df_res小，t-based q；方向验证通过。
- 未下载任何数据；mRNA从未升级为translation证据；未改docs/manifests/DECISIONS。
