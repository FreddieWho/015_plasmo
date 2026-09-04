# 疟原虫组成—功能项目 Agent 执行包

**项目代号：** PLASMODIUM-C2F  
**目标期刊：** Nature Microbiology（目标约束，不是必须维持的结论）  
**执行包版本：** v1.0  
**冻结日期：** 2026-09-02  
**当前里程碑：** M0｜立项与数据审计  
**当前主状态：** 翻译轴与可塑性轴仅并行观察；证据门 B 后必须二选一聚焦。

---

## 1. 本执行包的用途

本包不是一份“长提示词”，而是项目的持续状态与科学约束系统。其目标是让任意主 Agent 或子 Agent 在中断、换模型、换会话或并行施工后，仍能回答四个问题：

1. 项目的科学问题到底是什么；
2. 当前已经证明了什么、尚未证明什么；
3. 现在应完成哪个可验收交付物；
4. 哪些结果会导致继续、转向或停止。

**最高优先级原则：** 项目可以降级、转向或否证原假设，但不能靠增加算法、数据模态或湿实验来掩盖证据门失败。

## 2. Agent 必读顺序

每次新会话、断点恢复或更换主 Agent 时，按以下顺序读取：

1. `manifests/project_state.yaml`：当前唯一状态源；
2. `00_README_FIRST.md`：包结构与执行循环；
3. `01_AGENT_CHARTER_AND_NONDRIFT_CONTRACT.md`：不可静默修改的科学边界；
4. `02_SCIENTIFIC_SPECIFICATION.md`：科学假设、可证伪预测与主张等级；
5. `03_ROADMAP_GATES_AND_DELIVERABLES.md`：里程碑、证据门和交付物；
6. `04_ANALYSIS_CONTRACT_AND_MODEL_BOUNDARIES.md`：分析方法的最低要求与禁区；
7. 当前任务涉及的数据文档：`05_DATA_AUTO_DOWNLOAD.md` / `06_DATA_MANUAL_ACQUISITION.md`；
8. `manifests/gate_status.tsv`、`claim_registry.tsv`、`evidence_ledger.tsv`、最近的决策记录；
9. 其余专题文档按需读取。

## 3. 包内文档

| 文件 | 作用 |
|---|---|
| `01_AGENT_CHARTER_AND_NONDRIFT_CONTRACT.md` | 冻结科学边界、任务准入、变更控制和升级条件 |
| `02_SCIENTIFIC_SPECIFICATION.md` | 中心问题、H1–H5、主张梯子、允许与禁止的表述 |
| `03_ROADMAP_GATES_AND_DELIVERABLES.md` | WP0–WP4、M0–M5、Gate A–D、出口条件 |
| `04_ANALYSIS_CONTRACT_AND_MODEL_BOUNDARIES.md` | 系统发育、零模型、低复杂度、组成数据、模型复杂度边界 |
| `05_DATA_AUTO_DOWNLOAD.md` | Agent 可自动获取的数据、优先级、访问路径、下载审计规则 |
| `06_DATA_MANUAL_ACQUISITION.md` | 仅在自动获取失败时需要用户操作的数据；明确排除“联系作者”数据 |
| `07_LITERATURE_FRONTIER_AND_NOVELTY_MAP.md` | 高水平研究锚点、当前领域关切及避免重复的边界 |
| `08_EVIDENCE_CLAIM_AND_DECISION_PROTOCOL.md` | 证据账本、主张状态、Gate 评审与冲突处理 |
| `09_QA_REPRODUCIBILITY_AND_ACCEPTANCE.md` | 数据、统计、图件、代码和交付物的验收标准 |
| `10_MANUSCRIPT_STORYBOARD.md` | 目标论文的 Figure 1–5、备选故事和期刊降级标准 |
| `11_AGENT_MASTER_PROMPT.md` | 可直接交给主 Agent 的总控提示词 |
| `12_COLLABORATION_AND_EXPERT_INPUT.md` | 与寄生虫专家的输入接口、需要人工判断的节点 |

机器可读状态与模板位于 `manifests/` 和 `templates/`。原方案书与两篇起始论文位于 `sources/`。

## 4. 每次执行的固定循环

### 开始时

1. 读取 `project_state.yaml`；
2. 检查当前 milestone、active_gate、open_decisions、next_action；
3. 读取最近一次 session handoff 与 decision log；
4. 将拟执行任务映射到 **Aim → WP → Milestone → Gate → Deliverable**；
5. 若无法映射，任务不得开始。

### 执行时

- 原始数据只写入 `data/raw/`，派生结果只写入 `data/derived/`；
- 每个结论记录统计单位、对照/零模型、替代解释和敏感性结果；
- 每个新数据模态必须明确“它区分哪两个竞争解释”；
- 每个子 Agent 只能领取一个边界明确的任务包，并返回完整 provenance 与验收结果。

### 结束时

至少更新：

1. `project_state.yaml`；
2. `manifests/evidence_ledger.tsv`；
3. `manifests/claim_registry.tsv`；
4. `manifests/data_registry.tsv`（如发生数据变更）；
5. 一份 `templates/T03_SESSION_HANDOFF.md` 格式的交接记录；
6. 精确到文件和验收条件的下一步，不写“继续分析”。

## 5. 当前第一任务

主 Agent 启动后**不得直接开始全量下载或分析**。第一项交付物是：

> `M0_READINESS_AND_DATA_AUDIT.md`：验证执行包完整性，解析候选物种的当前参考组装，审计各数据源可访问性，给出发现层/深挖层/外群层的暂定面板及其风险，但不提前锁定主机制。

完成该审计并通过 M0 评审后，才允许进入 WP1。

## 6. 当前不允许做的事

- 从参考基因组 GC 含量直接宣称“突变偏倚”；
- 把基因数当成独立物种样本数；
- 用随机基因拆分得到的高 AUC 宣称跨物种泛化；
- 只用 RSCU、ENC、PCA 或 GO 富集支撑适应性；
- 在 Gate B 前大规模下载群体原始 reads、全量蛋白结构或宿主组学；
- 在 Gate D 前把湿实验写成既定任务；
- 同时长期维持翻译轴和可塑性轴两条等权主线；
- 为维持 Nature Microbiology 目标而隐瞒负结果或扩大叙事。

## 7. 成功的定义

成功不是“无论结果如何都投 Nature Microbiology”，而是完成一条可审计的推断链：

**真实组成现象 → 排除数据伪影 → 系统发育与序列背景反事实 → 稳定功能残差 → 独立功能证据 → 可证伪机制。**

若链条在某个 Gate 失败，按预设规则 pivot 或 stop，同样视为项目管理成功。
