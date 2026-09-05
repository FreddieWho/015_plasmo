# 证据、主张与决策协议

## 1. 三本账

### 数据账 `data_registry.tsv`
回答“数据从哪里来、版本是什么、能否支持当前问题”。

### 证据账 `evidence_ledger.tsv`
回答“具体结果是什么、质量如何、支持或削弱哪个主张”。

### 主张账 `claim_registry.tsv`
回答“当前允许说到什么程度、还缺什么证据”。

任何结论如果没有同时进入证据账和主张账，不得出现在汇报、摘要或主图。

## 2. 证据状态

| 状态 | 含义 |
|---|---|
| OBSERVED | 对已冻结数据的直接描述 |
| DERIVED | 通过可复现统计/计算得到 |
| INFERRED | 依赖模型假设的解释 |
| SUPPORTED | 由至少一类独立证据支持 |
| MECHANISTIC | 有决定性扰动或等价自然实验 |
| SPECULATIVE | 合理但未验证的假说 |
| REJECTED | 被数据或敏感性分析否证 |
| UNRESOLVED | 证据冲突或质量不足 |

## 3. 每条证据的必填字段

- evidence_id；
- 日期、Agent/分析版本；
- 输入数据与 checksum；
- 分析单元；
- 结果与效应量/不确定性；
- 零模型；
- 主要敏感性；
- 替代解释；
- 支持/削弱的 claim_id；
- Gate 影响；
- 可重复输出路径；
- 状态和限制。

## 4. Gate 评审规则

1. 由主 Agent 汇总，不允许单个子 Agent自行宣布 Gate 通过；
2. 同时展示支持证据、反证、质量缺口和最小下一步；
3. 先按预设标准判断，再讨论期刊叙事；
4. 决策只能是 GO、PIVOT、STOP/DOWNGRADE、HOLD；
5. HOLD 必须有一个成本受控且能改变决策的信息增益任务；
6. 用户和寄生虫专家对科学范围有最终决策权；
7. 决策写入 `gate_status.tsv` 和决策日志。

## 5. 冲突处理

### 数据冲突
优先检查 assembly/annotation、gene ID、阶段标签、批次、低复杂度与技术平台。不能通过平均相互矛盾结果来“达成一致”。

### 模型冲突
先比较各模型保持的约束是否不同；以预注册主模型为主，替代模型作为敏感性。不得选择最显著的模型作为最终模型。

### 文献冲突
区分物种、阶段、药物条件和测量层级。若确为矛盾，保留双方证据并设计可区分预测。

### 专家意见与计算冲突
专家意见用于识别错误假设和生物学边界，不直接覆盖数据。形成可检验的修订后再决策。

## 6. 负结果政策

- Gate 失败必须写入证据账；
- 不得把不显著改写成“趋势”而不报告区间；
- 不得在失败后无限细分亚组；
- 探索性新分组必须标记 POST_HOC，并在独立物种/数据中验证；
- 负结果可支持期刊降级、方法论文或可塑性转向。

## 7. 对外汇报最小格式

每个结果页必须包含：
1. 问题；
2. 统计单位；
3. 主结果与效应量；
4. 关键零模型/敏感性；
5. 允许的结论；
6. 不能得出的结论；
7. 对当前 Gate 的影响。

## 8. 会话交接

每次 session 必须用 `templates/T03_SESSION_HANDOFF.md`，至少写明：
- 已完成和未完成；
- 新增/变更文件；
- 失败尝试；
- 证据/主张更新；
- 当前硬阻塞；
- 下一任务及验收条件。

“下次继续分析”不合格。

## 9. M4R 证据规则（2026-09-05，D-036；不搞预注册）

- discovery 可自由探索，不设 outcome-blind 预注册门禁；独立验证优先（同一实验不同 modality = orthogonal evidence，不是 independent biological replication）。
- exploratory result 不因不是最优模型而自动删除；claim strength 与方法质量相匹配即可（dTE ratio = DISCOVERY/SUPPORTING，interaction 存活 = upgrade）。
- 防漂移规则（非预注册 Gate）：C 若 stage-corrected 后弱/阴且 ≥2 独立 ART 数据集无 Asn/LCR/regulator signal，则收为 boundary/negative result，不再加模态。
- 措辞禁令（违反即回退确认性约束）：禁止再写“deliberately designed bottleneck”“Asn is the first supply to collapse”“ART resistance is implemented through this axis”“acute mobilization/chronic conservation”作为已证机制；Asn supply 改用“Asn decoding demand is unusually high / Asn-related translational supply may constitute a constraint / starvation-associated tRNA data are consistent-inconsistent with this hypothesis”（由数据定）；dTE ratio 保留探索标签；acute/chronic 反号只作 hypothesis。
