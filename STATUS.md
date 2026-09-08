# STATUS — 项目唯一入口（两部分同时更新，不矛盾）

## 给人读的进展

已经完成的是路线重启。旧的跨物种密码子故事在去年验证失败后正式封存、永久保留；新发现（调控蛋白富含 Asn、组成驱动的蛋白变化与耐药位点走的是两条路）让用户批准了一次有限重启：新问题是基因组组成变化如何传到蛋白序列、落在蛋白的哪些部位、什么时候被生命周期或药物应激利用起来。

正在做的是手稿打包收尾。M5 终审已出（D-039）：D 进化主骨架 CONDITIONAL-STRONG（LOO 7/7、D2b 框内存活、D6 规则满足、Pf8 群体边界）；A 富集臂正式集复现但 state 臂失败（liver 检测推翻、zygote null、不收候选），Figure 4 有界使用、不宣称 co-option；C 为 Figure 5/Extended 第二语境。T01 不激活。L-010 已落盘（curated 欠定，T01 不触发）；Fig1–5 共 14 panels 多模态目检通过（含 8 项修复 + Fig5a 诚实性修复）；主文本草稿 v0.1 + 主表 S1–S6 已组装。

卡在三处。祖先序列重建对部分位点不可靠，只能标记后排除；急性药物×基因型的翻译层公共数据根本不存在，这是已声明的缺口、不用弱替代品硬填；重复功能新颖性已被 2026 年一篇 Nature Communications 抢先，标题必须绕开它。

准备这样解决。M5 终审（D-039）已裁决 D+A+C 三线各安其位：D 主骨架 Fig2–3、A 有界 Fig4（不宣称 co-option）、C 为 Fig5/Extended 第二语境，T01 不激活。下一步是手稿打包；投稿前必做只有 L-010 真 GO 复核，SNP-level burden 已明确 deferred，不再开新分析分支。

```
2026/9/5
ROADMAP  [########] 8/8 节点（R0–R5 完成；M5 终审 D+A+C 各安其位，T01 不激活）
本周投入  科学问题 ████████░░ 80%   基础设施 ██░░░░░░░░ 20%

偏离程度  低
偏离位置  M4R 是用户批准的 bounded rescue（D-036），旧 Gate C PIVOT 保留有效；偏离的不是旧假说而是新问题 composition→protein→function。
建议     冻结新分析分支（C 维持、Puf1/GAP-4 维持缺口声明）；只做手稿打包 + L-010 复核。
```

## 给 agent 的接手信息

- 活跃节点：M4R/WP4R（M4R_RESCUE_ACTIVE，backbone COMPOSITION_TO_PROTEIN_FUNCTION，routes D+A+C）
- 核心文件：docs/M4R_ROUTE_RESET_20260905.md、docs/tasks/M4R-{D,A,C}_TASK_PACKET.md、docs/manifests/project_state.yaml
- 复现：sha256sum -c data/checksums/*.sha256；L2 v2 见 data/derived/WP4/L2_site_composition/；M4R 输出 data/derived/WP4R/M4R-{D,A,C}/
- 最近 DECISIONS：D-039 M5终审（D+A+C 各安其位；偏离预注册字面但服从立法本意的理由见终审 §2）
- 下一步：手稿已打包（manuscript/：14 panels + legends + 主表 S1–S6 + 主文本 v0.1）；投稿前剩余：用户定稿标题/摘要 + 期刊格式适配
