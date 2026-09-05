# M4R Round-1 对抗性审阅（2026-09-05，主Agent自审）

**审阅对象：** M4R-0 package reset + 三线 round-1 结果 + M5 convergence。
**审阅标准：** 是否偏离设定科学目标（composition→protein→function）；结论与数据一致性；账本完整性；需要补充/修复项。

## 一、与科学目标的偏移评估：无实质偏移（低）

- 新问题（composition→protein→function）未被偷偷换回旧 codon programme；Gate C PIVOT 与全部负结果原样保留；三线输出均停在 association，未升 mechanism；措辞禁令（supply collapse / dual-axes / mobilization-conservation）在全部报告与账本中未见违反。
- T01 降级未被删除，标记 SUPERSEDED_FOR_NOW 并保留重激活条件——合规。
- 唯一叙事张力：M5 报告说 "D+A leading"。审阅确认这是证据排序（D backbone 强 + A 富集强），不是为期刊选答案；报告同时声明 L5/L6 未补齐。可接受。

## 二、进度核对：round-1 六项交付齐（§20）

ROADMAP 已更新为 8 节点（R0–R2 完成、R3 阴性、R4R round-1 完成、R5 待收敛）；STATUS 进度条已同步修正（原 5/10 与 ROADMAP 不对应，属格式杜撰，已改 6/8）。PLAN.md 未动（合规，agent 无权修改）。

## 三、科学性弱点（按严重度排序）

1. **D2 "99% structured churn" 部分是同义反复（最重要）。** core-181 的构造就排除了 LCR 重/快进化蛋白（D6 自己发现 AP2/PUF 为 0 个入选）——在这个采样框里说"churn 不在 LCR 里"是被构造出来的。且项目约定 LCR caller 在蛋白上 0%/0%（不可用）。M5 报告写 "structured churn" 的强度需要下调为："在保守单拷贝核心基因内，耦合 turnover 不依赖 repeat 序列"；并补一条全蛋白组（不限 core-OG）的 LCR/IDR 归属作为 round-2 必做项，否则 Figure 3 的 architecture 论断站不住。
2. **C 的 stage-adjusted 结论强于其实际效力。** stage reference 只解释 ~1e-4–3% 方差；"stage 不是驱动者"只在"这个弱参考下信号不变"意义上成立，不能排除更强 stage 定义下的混杂。且 STAGEADJ 后 27 组中 6 组原显著转弱/转负（Cam3II 24/40h、revWT 8/40/48h、Dd2_R539T 48h），报告"存活"是选择性陈述——应改为"16/27 存活、Dd2 背景稳定、Cam3II 混合"。同时 stagevec 单一 ring-schizont 轴是粗糙代理。C 的 CONDITIONAL 判定不变，但措辞必须收紧。
3. **D3 matched controls 未达 packet 要求。** 只匹配了 conservation；packet 要求再匹配 protein/domain、structural context、essentiality、alignment quality。弱版措辞保住了结论，但"matched"字样在 M5 中给读者过高预期。round-2 补 domain/essentiality 维度的匹配或改名 "conservation-matched only"。
4. **PUF 口径不一致。** D5 的 PUF 信号（OR=3.03）用的是 56 基因 PUF_RNA 宽集；A 的严格 PUF 集只有 2 基因不可判。CLM09 与 M5 报告里的 "AP2/PUF" 应改为 "AP2（强）+ 宽 RNA-binding（弱）；严格 PUF 未定"。
5. **A 的正式 enrichment 实为 product-description regex 衍生**（GFF 无 GO 条目），不是 packet 要求的 GO/InterPro 正式集。"正式"字样过头；round-2 的 InterPro 复现是硬项，不是可选项。
6. **C 的 acute 三个"独立"数据集**：scRNA（isogenic K13）与 GSE151189 真独立，但都是 mRNA；persistence GSE225340 是 dormancy 而非 acute DHA。cross-dataset independence 维度的评分在 M5 中偏高半档。
7. **小项：** D6 pilot 的 4/4 方向一致无统计意义（n=4 的"方向一致"不构成证据，报告已标 UNRESOLVED，保持）；dTE interaction 的 t 检验 df 过小（规则 2 处理正确）；GSE59099 halflife 双键合并方法合理但需写进 claim_impact（已写 NOTE 被删除——应在 claim_impact 补一句）。

## 四、工程/账本缺口（本轮已修复项标注 ✓）

- ✓ **git 未提交**（整个 M4R 包+manifests 全脏）→ 已提交 b456e3f。data/derived 按 D-034 惯例不入 git（sha256 溯源在 input_manifest）。
- ✓ **PACKAGE_MANIFEST.sha256 59 项失效** → 已重生成，0 项失配。
- ✓ **infra/bioinf-data-index 缺失**（全局约束：新增外部生信数据交付前必须更新）→ 已建 `infra/bioinf-data-index/M4R_20260905.md`（9 个新数据集索引+汇总）。
- ✓ **docs/CHANGELOG.md 停在 v1.0** → 已补 v1.1-M4R。
- ✓ **ROADMAP.md 无 M4R 节点** → 已补 R4R；R3 标阴性完成。
- ✗ **TODO.md 残留过时项**（"M3 供需 R3"、"T01 期刊降级签字"、"M4 Pf 机制聚焦"等 Gate C PIVOT 前的勾选项未清理；缺 AGENTS.md 要求的"分支记录"段）→ 需清理但保留变更记录。
- ✗ **gate_status.tsv 用了非标准状态值** OPEN_M4R_RESCUE_ACTIVE（08 §4 只允许 GO/PIVOT/STOP/HOLD）→ 改为 HOLD（HOLD 的信息增益任务 = M4R round-2）或扩协议定义；建议改 HOLD。
- ✗ **LEADS.md 未吸收 round-1 产生的支线**（如：D2 全蛋白组归属、PUF 严格集构建、branch-LOO、dTE interaction 正式 anota2seq）→ 建议补 3–4 条。

## 五、需要补充的分析（round-2 硬项，不扩主题）

1. **D2b 全蛋白组归属**（不限 core-181）：coupled-turnover 候选的全蛋白 LCR/IDR/domain 归属 + 采样框偏差声明。这是 Figure 3 的生死项。
2. **D6 形式化**：≥2 独立 transitions 的预注册 pole 对比（不是 mafft pilot n=4）。
3. **A 的 InterPro 正式富集 + 严格 PUF 集构建**（从 PlasmoDB annotation/orthology 而不是 regex）。
4. **C 措辞修正**：stage-adjusted 弱参考声明 + 16/27 存活如实表述 + halflife 合并方法入 claim_impact。
5. branch-LOO（D1 的系统发育敏感性）与 Pf8 罕见变异边界可缓，但在投稿前必须。

## 六、对 M5 判定的复核

D+A leading 的排序本身成立，但图级判断应降为：**D backbone 可立（Figure 2–3），前提是 D2b 全蛋白组归属在同向存活**；若 D2b 显示 churn 实为 LCR 富集（采样框外），Figure 3 与 general-principle 强度都要降。A/C 维持 round-1 判定。T01 不激活，但 round-2 后必须重估。

**结论：无目标偏移，无账本断链；1 个实质科学弱点（D2 采样框）+ 2 个措辞过强（C stage、A"正式"）需要修正后 round-2 才可启动。**
