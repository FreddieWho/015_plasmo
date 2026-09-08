# TODO — 当前颗粒度待办（按执行顺序，重要的在前）

- [x] M4R-0 状态与证据修复 + package/state reset（2026-09-05，D-036：01/02/03/04/07/08/10/11 修订；state 1.1-M4R；CLM04 修订+CLM08/09/10 新增；gate M4R_CONVERGENCE；EVID-M4R-000；L33-35；ROUTE_RESET 权威文档；T01 标记 SUPERSEDED_FOR_NOW）
- [x] M4R-D/A/C 第一轮信息增益分析（2026-09-05：D/A worker回包均为CONDITIONAL；C两轮超时后主Agent收尾完成，CONDITIONAL；MWU exact-hang根因修复）
    - [x] D：branch-aware L2 + accessibility + architecture + matched controls + regulatory enrichment + 自然实验
    - [x] A：正式 enrichment + lifecycle 映射 + 候选收缩（本轮不收shortlist）
    - [x] C：stage-adjusted DHA + acute/chronic + feature 分解 + 独立 ART 验证 + dTE upgrade
- [x] M5 convergence review round-1（2026-09-05，D-037：D+A leading，D+C第二context，降级暂不执行；见docs/M5_CONVERGENCE_REPORT_M4R_ROUND1.md）
- [x] M4R round-2 bounded增量（2026-09-05，D-038启动，双包并行一次通过无超时）
    - [x] D2b 全蛋白组归属（SURVIVES框内；harbor匹配后LCR零超额；采样框偏差定量）
    - [x] D6形式化（T2+T3规则满足，T1注册阴性，CHROM null）+ branch-LOO 7/7 + D3多维匹配（改名）
    - [x] Pf8群体边界（34/34 markers分离+CNV率；SNP-level deferred）
    - [x] A lifecycle映射（三极：gametocyte弱/zygote null/liver检测推翻排除）+ PfAP2-P/PbApiAP2表达 + InterPro正式集 + 严格PUF集（n=2不可判）；不收shortlist
- [x] M5终审（2026-09-05，D-039：D+A+C三线各安其位，T01不激活；见docs/M5_FINAL_REVIEW_20260905.md）
- [x] 手稿打包（2026-09-08 用户指令）：L-010 真 GO 复核落盘（curated 欠定，T01 不触发）+ Fig1–5 图包 14 panels 多模态目检（8 项修复 + Fig5a dir-sig 诚实性修复，全部复检通过）+ 主文本草稿 v0.1 + 主表 S1–S6 组装

- [x] M4R T1 数据准备（2026-09-05：GEO 六集 + GCN5 全套 VERIFIED，~540MB；essentiality/K13scRNA 转人工）
    - [x] GEO：GSE75795/GSE120448/GSE134268/GSE120488/GSE225340/GSE59099 落盘+manifest+sha256
    - [x] GCN5 2026（PMC13438671，MOESM1-7，Source Data 在手；DOI 已验）
    - [x] 人工：essentiality 三篇 30 文件 + K13 scRNA zip 全部到货验收入库（2026-09-05，inbox+raw+manifest+sha256，registry 已翻 VERIFIED）
- [-] M3 供需 R3（D-008 分层）   下载子项可并行，需求计算阻塞于下载完成  【已废止：Gate C PIVOT 2026-09-04 / M4R reset 2026-09-05，见分支记录】
    - [x] 主：D006/D011/D012/D013 落盘登记 VERIFIED（2026-09-04，21MB+1.3MB）
    - [-] 辅：D008 M03 部分到货（正文+supp PDF 已入库，xlsx 待定，不阻塞）+ D014/Pv4 按需
    - [-] D_{c,s,t} 需求计算，主看 Pf-vs-Pk 超越 genome_gc，辅线只做支持
    - [-] tRNA/ribo/mRNA 第二正交层整合与跨谱系一致性
- [x] D_{c,s,t} 主判定跑通（2026-09-04，8 TSVs + 交付段全齐）
- [x] M3-03 跑完判阴（2026-09-04，0 early hits，assay-sensitive）→ Gate C PIVOT
- [x] 攻击性审查修复 1-4（2026-09-04，D-033：git+环境锁 / GPL18893 溯源 / tRNA 隔离 / I5 双轨 0=0 spearman=0.9991 PIVOT 稳健）
    - [x] git init + 环境锁文件（提交 015f422/54ee48d/1590a42；LFS 迁移前 fbd2770/37944a8/aa64c46）
    - [x] M3-03 溯源修复 + GPL18893 入 data_registry
    - [x] 隔离 tRNA supply/charging 表（NOT_FOR_CLAIM）
    - [x] I5 双轨补跑（遮蔽轨命中 0=0，判定稳健）
- [~] L1/L2 双方向（D-018/D-019/D-020，方案已交付，L1 后台运行中 b8a06ba14）
    - [x] 文献与方法核查（2026-09-04：L1 空位确认；L2 组合角度新）
    - [x] L1 方案 + L2 方案（L1_PLAN_CODONxSTRESS.md / L2_PLAN_SITE_COMPOSITION.md）
    - [x] L1 结果验收（2026-09-04 全扫版 DOUBLE_POSITIVE，D-023：AAT/Val/Ala 三家族双阳，Lys/Ile 未复现）
    - [x] L2 v1 实施+验收（2026-09-04，D-026：H-L2a 反向、H-L2b 成立；官方 orthogroup 到货后复跑）
- [x] L1b 解释（2026-09-04，D-025：5→1，AAT 唯一超组成双阳，Val/Ala=第一位回声）
- [-] T01 期刊降级签字（B 档三选一）        ← 阻塞 M4 收口口径  【已废止：Gate C PIVOT 2026-09-04 / M4R reset 2026-09-05，见分支记录】
- [-] M4 Pf 机制聚焦（AAA/s2U/K13 + 2 竞争模型 + 逃逸边界，在手数据）  【已废止：Gate C PIVOT 2026-09-04 / M4R reset 2026-09-05，见分支记录】
- [x] M3-02 功能层验收（2026-09-04，9 文件 + GSE151189）
- [-] M3-01 + M3-02 会合：LysAAA/AAG + s2U 轴第二层检验 → Gate C 证据卡  【已废止：Gate C PIVOT 2026-09-04 / M4R reset 2026-09-05，见分支记录】
- [-] 定 Pv P01 vs Salvador-I 双跑敏感性方案  【已废止：Gate C PIVOT 2026-09-04 / M4R reset 2026-09-05，见分支记录】
- [-] Gate C 评审材料组装（每个候选跨层证据卡 + 冲突标记 UNRESOLVED）  【已废止：Gate C PIVOT 2026-09-04 / M4R reset 2026-09-05，见分支记录】

## 变更记录

- 2026-09-05：人工包入库完成（data/tmp 已清空：essentiality 30 xlsx 按三篇分目录进 manual_inbox；K13 zip 213MB 进 inbox + 55 文件解压至 data/raw/k13scrna/，14 DGE GENE×cell 已验 PF3D7；Elsworth Data_S3=核心 essentiality 表；manifest+sha256 77 项全 OK；registry M4R-D02/C03 翻 VERIFIED）。
- 2026-09-05：M4R T1 数据准备完成（M4R-A02/A03×3/C01/C02b/A01 共7集 VERIFIED，manifest+sha256+m4r-t1；data_registry 追加9行；ChIP trio 纠正为单研究 PMID 32198457；essentiality 因 PMC PoW/Cloudflare 转人工 M4R-T3-04；K13 scRNA 验号 PRJNA1049964+Zenodo 20344254 待取）。原因：M4R 搜索包落地后的数据准备延续。

- 2026-09-04：M3 系列完结、Gate C PIVOT（D-016/D-017）；新增攻击性审查修复四件、T01 签字、L1/L2 双方向（D-018）。原因：M3-03 阴性触发自动 PIVOT。
- 2026-09-04：L1 实施并验收（DISCOVERY_ONLY，D-021）；D-019 确立探索性研究不搞预注册；D-022 全空间扫描纪律（入 ~/AGENTS.md），L1 复现层扩为 61 密码子全扫重跑中；D-020 T01 延至 L1/L2 结果齐后。原因：用户三次指示。
- [x] 环3 蛋白层检验（D-029 D008 未确认 + D-031 Li 饥饿蛋白组部分阳性 rho=+0.167）
- [x] Mok 2021 MOESM6 基因型情境检验（2026-09-04，D-032：反号——抗性稳态节约 Asn 富集蛋白，与急性饥饿动员互补，D-028 双轴论点获第三语境支持）
- [x] L1f 翻译层第三数据集验证（2026-09-04，D-030：GSE58402 Ribo-seq，生长期条件性支持，19h rho=+0.35 p~1e-24）
- [x] L2 硬化 v2（2026-09-04，D-028：Grantham+QC+ML 树全部落地，v1 结论存活；官方 OG 仍缺，已披露为小局限）
- 2026-09-04：攻击性审查修复 1-4 全部完成（D-033），四子项逐个打勾。原因：用户指示"完成修复"；I5 双轨补跑证实 PIVOT 判定对 LCR 遮蔽稳健（0=0，rho=0.9991）。

## 变更记录（追加，只加不删）

- 2026-09-05：M4R rescue 启动（D-036 用户授权）：M4R-0 完成 + package 1.1-M4R + state（M4R/WP4R/COMPOSITION_TO_PROTEIN_FUNCTION/D+A+C/M4R_RESCUE_ACTIVE）+ claim 重审计（CLM08/09/10）+ ROUTE_RESET 权威文档 + M4R-D/A/C 三包并行。原因：M4 Asn/L2 新发现 + search supplement 三线锚点；旧 Gate C PIVOT 保留不推翻。
- 2026-09-05：M4R round-1完成（D-037）：D/A worker回包CONDITIONAL + C主Agent收尾CONDITIONAL（修MWU exact-hang）；M5 convergence判D+A leading；EVID-M4R-D/A/C-001入账。原因：三线第一轮是M4R→M5的约定收敛点；C两轮worker超时须如实记录并换主Agent有界收尾。
- 2026-09-05：对抗性审阅完成（docs/M4R_ADVERSARIAL_REVIEW_20260905.md）：无目标偏移；修复 git 提交/PACKAGE_MANIFEST 失效 59 项/bioinf-data-index 缺失/CHANGELOG v1.1/ROADMAP R4R/STATUS 进度条；gate M4R_CONVERGENCE 状态改 HOLD（信息增益任务=round-2）。待办：D2b 全蛋白组归属（Fig3 生死项）、D6 形式化、A InterPro 正式集、C 措辞修正（16/27 存活+弱 stage 参考）。

## 分支记录

- **当前执行分支：** M4R/WP4R（D backbone + A/C 并行 → M5 convergence），D-036/D-037。
- **停止分支：** 跨谱系 codon-stress programme（Gate C PIVOT 2026-09-04，永久）；Pf AAA/s2U/K13 唯一 M4 主线（2026-09-05 降 boundary context）；上列 6 条 [-] 已废止项随分支停止。
- **暂缓分支：** T01 期刊降级（SUPERSEDED_FOR_NOW，M5 终审时重估）；Pv 双株敏感性（MOOT）。
- **并行分支：** M4R-D / M4R-A / M4R-C 三线（round-1 完成，round-2 bounded）。
- 2026-09-05：M4R round-2完成+M5终审（D-038/D-039）：D/A双包一次通过无超时；D CONDITIONAL-STRONG（LOO 7/7、D2b框内SURVIVES、D6规则满足、Pf8边界），A CONDITIONAL不变（AP2正式集复现35.3、liver排除、zygote null、无shortlist）；终审D+A+C各安其位、T01不激活（偏离预注册字面理由见终审§2）；EVID-M4R-D/A-002入账；LEADS L-006/L-008并入、L-007部分、L-010新增。
- 2026-09-08：手稿打包完成（D-041）：L-010 落盘（curated 欠定，T01 不触发）+ Fig1–5 14 panels 目检（8 项修复 + Fig5a dir-sig 诚实性修复，全部复检通过）+ Fig4a GO 行并入 + 主文本 v0.1 + 主表 S1–S6 + FIG_FIXLOG 归档。原因：用户指令打包 + 打包前多模态核验。
