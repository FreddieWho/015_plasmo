# M4R-C ROUND1 outputs (2026-09-05)
问题：stage/confounder修正后ART signal剩多少？acute vs chronic是否独立成立？
方法：spearman + length-adjusted + Asn-top10 MW + noSchiz双轨；BH per contrast。
主表：M4RC_feature_decomposition.tsv（142 contrasts×10 features=1420行）。
分区：M4RC_acute_chronic_partition.tsv（asn_frac行+arm标签）。
dTE：M4RC_dTE_interaction.tsv + M4RC_dTE_upgrade.json（规则2：同向弱支持）。
脚本：run_m4rc.py（worker初版）/ run_m4rc_salvage.py（STEP版）/ run_m4rc_mine.py（主Agent收尾版，MWU asymptotic修复）。
输入版本见 input_manifest.tsv（本目录待补：见下）。
报告：M4R-C_ROUND1_REPORT.md（CONDITIONAL）+ claim_impact.md。
