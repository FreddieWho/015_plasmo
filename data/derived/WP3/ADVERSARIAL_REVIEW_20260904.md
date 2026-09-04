# 攻击性审查 — M3 实现与结果（2026-09-04）

审查范围：M3-01 demand、M3-02 功能提取、M3-03 DHA×K13、Gate C PIVOT 决策本身。
方法：逐条实证核查（文件/代码/输出已验证），不按汇报口径。

---

## A. 结果/决策层（最重）

### A1. M3-03 的阴性证据权重被高估 —「assay-sensitive null」说法不成立
- 需求份额是**组成型统计量**（各行和=1），被高表达看家基因主导。双色芯片经 lowess 归一化（前提：大多数基因不变），**归一化本身就把需求份额的差值压向零**。即使存在真实的特定转录本翻译重编程，聚合份额也几乎不动 → 0 hits 在设计上大概率必然出现。
- 所谓阳性对照（K13 mRNA +0.4..+1.0 log2）检验的是**另一个统计量**（单基因 mRNA），不能证明需求份额统计量有检测力；且晚时间点 K13 上调很可能就是**发育延迟混杂**（DHA 拖慢周期 → 阶段偏移 → K13 表观上调），早期 8h 只有 n=1。
- 结论：M3-03 应改述为「**统计量按构造不敏感 + 功效仅够大效应（3v3）的阴性**」，不是「灵敏检测下的阴性」。GATE_C_REVIEW.md 的 assay-sensitive 措辞需修订。
- 连锁影响：PIVOT 本身不倒——真正的 GO 否决票是 M3-01 预注册主判定 13/36（同样低功效：符号检验 + 需求差 1e-2 量级），HOLD 任务已完成流程义务。但「两线失败」里 M3-03 那一线的权重必须下调。

### A2. CLM04 的 Pf 域支撑有循环风险：M3-02 是文献再提取，不是本项目新证据
- M3-02 的核心事实（s2U↓、LysAAA 驱动、K13 52/57 AAA、PfMnmA cKD 4%→11%）全部出自 Small-Saunders 2024 本人数据；M3-01 的 TE/decay 关联与 Li 2024（anota2seq polysome/total、actD decay、密码子层）**方法同构**，属验证性复现。
- 即：CLM04「Pf 内 ≥2 正交层」里，一层是别人的论文，一层是别人论文的复现。本项目独有贡献 = 组成定量框架 + 阴性边界。
- CLM04 的 allowed_language 必须收紧为「与已发表机制一致、且在组成框架下定量化」，禁止暗示本项目独立发现了 Pf 翻译应激机制。否则审稿人一句「你在复述 Small-Saunders」即致命（07 §3 已预警）。

### A3. 跨物种 13/36 是「无证据」不是「证伪」
- 已部分披露（D-013），但 STATUS/GATE_C 行文中「失败」字样偏强。符号检验 + 亚家族拆分后功效进一步稀释。正确表述：跨谱系程序**未获支持**，不排除低功效漏检。这不改变 PIVOT（GO 举证责任在主张方），但影响摘要措辞。

## B. 实现层

### B1. I5（LCR 双轨）在 M3-03 完全缺席，M3-01 半缺席 — 冻结不变量违约
- `run_m3_03_dha.py` 无任何 masked/codm_ 轨道（grep 实证：0 处）。
- M3-01 的 demand 表与 expression-weighted null 用未屏蔽计数（`pre="cod_"`），却与 M2 的 **masked** 指标（`prop_explained_GC3_masked`）直接做差比较 —— 双轨基线不一致。
- 预期影响小（M2 双轨差 ~1e-3），但这是冻结不变量，要么补跑要么记正式 deviation。

### B2. M3-03 溯源断链：input_manifest 指向 /tmp/m3plat.tsv
- `/tmp` 是易失目录，manifest 记录的 sha 对应文件不在项目内。应指向 `data/raw/geo/GPL18893/GPL18893_family.soft.gz` + 提取步骤。
- 且 GPL18893 soft（259MB）进了 raw/ 却无 manifest、无 sha、未入 data_registry —— 违反 05 §1 下载合同。

### B3. tRNA 供应层建在人类（HS_）命名空间上，且 charging 层静默缺失
- GEO 计数文件实证为 `HS_Ala-AGC-1-1.1` 命名（人源 tRNA 参考）。Pf 仅 45 个 isoacceptor，用人类 tRNA 标签聚合出的 S/D 表生物学有效性未知。
- charging 正则失配 → `M3-01_tRNA_charging.tsv` 空文件，仅 1 字节，无显式报错。
- 处置：供应/charging 两表应隔离标记 NOT_FOR_CLAIM，Gate C/D 证据卡不得引用。

### B4. 全项目无 git — 09 §1 硬性违规
- `git rev-parse` 实证 not a git repository。所有派生产物、脚本、账本零版本控制。最终复现包（09 §9）还要求环境锁文件，目前只有 D-011 的 LD_PRELOAD 土办法。

### B5. TE/decay 基因级 p 值未做基因聚类/bootstrap
- 04 §8 要求基因集层面按 orthogroup/物种 bootstrap；M3-01 把 5168 基因当独立点（n=5168, p 至 1e-108）。密码子偏好由氨基酸组成和基因长度强结构化，p 值系统性膨胀。方向性结论大概率稳健，但「271/405 显著」的计数不可按面值引用。

### B6. M3-03 输出表卫生问题
- Cam3II n=1 分层产生 `p=0.0e+0, q=NA` 行混入交付 TSV（实证在 M3-03_DHA_effect.tsv / K13_DiD.tsv）；DiD 用 z 近似 + 3v3 var(ddof=1) 近零，无置换退路。负结果合规要求「报告区间」，这些行目前是不可解释的占位数值。

## C. 站得住的部分（攻击后存活）

- M1/M2 主体：梯度、84.5%/78.7% 定量、LOO、双轨（在 M2 层是做了的）、伪影筛查、文献边界（D-014/015）——这是全项目最硬、最新颖的部分。
- 数据治理：manifest/checksum/processed-first/manual_inbox 协议全程合规；负结果入账未删；D-016 按预注册原文执行。
- PIVOT 的**程序正义**成立：GO 举证未达标 + HOLD 任务已交付 + 预注册规则触发自动落 B。

## D. 修复优先级（已入 TODO）

1. [必须/廉价] git init + 环境锁（把 D-011 落成 environment 文件）
2. [必须/廉价] M3-03 溯源修复 + GPL18893 入 data_registry
3. [必须] 隔离 tRNA supply/charging 表（NOT_FOR_CLAIM）
4. [应当] I5 双轨补跑（M3-01 demand/M3-03）或正式 deviation 记录
5. [应当] GATE_C_REVIEW 措辞修订（A1/A2 两条）+ CLM04 allowed_language 收紧
6. [用户定] M4 改名「边界与整合」还是保留「机制聚焦」
7. [不建议现在做] 基因级 DHA×AAA 交互才是功效正确的检验，但属 POST_HOC；若做必须独立验证。当前不动。

——审查人：主 agent 自查（2026-09-04），实证核查见各条引用路径
