# ROADMAP — PLASMODIUM-C2F 技术路径（每个节点绑定 PLAN 假设）

> ROADMAP 宏观稳定；微观高频变动放 TODO.md，同一件事不两处重复。

- [x] **R0 M0/WP0 立项与数据审计** → 检验：能否提出 H1（数据是否够格）。状态：完成，G0 GO 2026-09-03。含 `[infra]` 子工作：D001–D005 落盘与 checksum。
- [x] **R1 M1/WP1 组成演化地图** → 检验：H1（多状态转换是否存在）。状态：完成，Gate A GO。GC 18.2%→52.3%，独立性 Grade A，7/7 LOO。
- [x] **R2 M2/WP2 跨层反事实分解** → 检验：H2（机械传导 vs 残差）。状态：完成，Gate B GO（窄化）。密码子 84.5% 被 GC3 解释、氨基酸 CLR 78.7% 被 genome_gc 解释；7 solid + L provisional，13 GC-coupled。
- [ ] **R3 M3/WP3 功能与扰动整合（D-008 分层：主 Pf-vs-Pk，辅 Pb/Pv）** → 检验：H3（残差是否预测翻译供需/扰动）。输入 D006–D013 processed-first，做 D_{c,s,t} 表达加权需求 + tRNA/ribo/mRNA 第二层，主判定只看 Pf-vs-Pk。Gate C：至少一模块 ≥2 正交证据。
- [ ] **R4 M4/WP4 机制聚焦与转化边界** → 检验：H4（条件脆弱性）/ H5 备选。收缩到 1 主机制 + 1–3 节点，明确两个竞争机制，定 Gate D（不做 / 最小实验 / 降级）。
- [ ] **R5 M5 收口投稿** → 按证据定 A 目标故事 / B 资源桥接 / C 中性传导 / D 可塑性转向。不为期刊追加模态。

备注：纯基础设施节点标 `[infra]`；若连续出现两个 `[infra]` 节点必须提醒用户。本路段暂无连续 infra。
