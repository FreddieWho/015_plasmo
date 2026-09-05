# ROADMAP — PLASMODIUM-C2F 技术路径（每个节点绑定 PLAN 假设）

> ROADMAP 宏观稳定；微观高频变动放 TODO.md，同一件事不两处重复。

- [x] **R0 M0/WP0 立项与数据审计** → 检验：能否提出 H1（数据是否够格）。状态：完成，G0 GO 2026-09-03。含 `[infra]` 子工作：D001–D005 落盘与 checksum。
- [x] **R1 M1/WP1 组成演化地图** → 检验：H1（多状态转换是否存在）。状态：完成，Gate A GO。GC 18.2%→52.3%，独立性 Grade A，7/7 LOO。
- [x] **R2 M2/WP2 跨层反事实分解** → 检验：H2（机械传导 vs 残差）。状态：完成，Gate B GO（窄化）。密码子 84.5% 被 GC3 解释、氨基酸 CLR 78.7% 被 genome_gc 解释；7 solid + L provisional，13 GC-coupled。
- [ ] **R3 M3/WP3 功能与扰动整合（D-008 分层：主 Pf-vs-Pk，辅 Pb/Pv）** → 检验：H3（残差是否预测翻译供需/扰动）。输入 D006–D013 processed-first，做 D_{c,s,t} 表达加权需求 + tRNA/ribo/mRNA 第二层，主判定只看 Pf-vs-Pk。Gate C：至少一模块 ≥2 正交证据。状态：完成（阴），Gate C PIVOT 2026-09-04——H3 跨谱系版死亡，永久保留为负结果。
- [x] **R4R M4R/WP4R bounded rescue（D-036）** → round-1（D-037）+ round-2（D-038）完成：D CONDITIONAL-STRONG / A CONDITIONAL / C CONDITIONAL。
- [x] **R5 M5 收敛（D-039，2026-09-05）** → D+A+C 三线各安其位（D backbone Fig2–3；A 有界 Fig4；C Fig5/Extended）；T01 不激活。下一步：手稿打包（package figures + 方法文本），L-008/L-010 为投稿前必做/待挖掘。

备注：纯基础设施节点标 `[infra]`；若连续出现两个 `[infra]` 节点必须提醒用户。本路段暂无连续 infra。
