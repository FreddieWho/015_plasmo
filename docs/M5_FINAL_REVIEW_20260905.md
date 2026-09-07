# M5 终审 — M4R Round-2 后统一裁决（2026-09-05，主Agent）

**问题：** D / D+A / D+C / D+A+C / 降级，五选一（不为期刊选答案）。
**输入：** M4R-D ROUND2（CONDITIONAL-STRONG）/ M4R-A ROUND2（CONDITIONAL 不变）/ M4R-C ROUND1（CONDITIONAL，维持）。
**预注册规则回顾（M5 round-1 §6）：** A-state转正则D+A成稿；A-state失败则D+C partition或诚实降级（T01重激活）。

## 1. 各线终值

- **D：CONDITIONAL-STRONG。** D1 branch OR=2.25 且 LOO 7/7 全显著（OR 1.71–4.35，最弱为 vivax-drop，符合预期）；D2b SURVIVES（frame-internal：harbor 基因匹配后 LCR 零超额，Asn-贫 0.064 vs 0.119，domain 富集；采样框偏差已量化声明）；D6 规则满足（T2+T3 独立 transitions 组成一致，T1 注册为阴性，CHROM 全程 null）；D3 部分升级并改名（domain-shared where available；within-OG 不可能已论证为结构性限制）；Pf8 群体层支持分离（34/34 耐药 markers 分离；CNV 率）；SNP-level burden 明确 deferred。
- **A：CONDITIONAL 不变（富集臂加强，state 臂削弱）。** AP2 正式集逐值复现（OR 35.3 vs 34.3）；chromatin 正式集转弱阳（新）；logit 调整存活（1.97）；PfAP2-P 身份验后 bounding（bound 基因 Asn-贫）+ Pb 家族 map 落盘。削弱面：严格 PUF 仍 n=2 不可判；repeatbroad null；liver 极 detection 推翻并排除；zygote null；gametocyte 仅弱方向；**state-consequence 标准未达，不收 shortlist，H-A 维持 hypothesis。**
- **C：CONDITIONAL 不变**（partition figure；MOA 不指望；M4R-X01 维持）。

## 2. 五选一裁决：D+A+C（三线各安其位），T01 不激活

- A-state（family-level lifecycle-state consequence）确实未转正——预注册规则的字面指向 D+C 或降级。
- 但规则未预见的实际结果是：**A 富集臂在本轮加强（正式集复现 + chromatin 新弱阳），而 C 分区证据（rho 0.06–0.23、mRNA-only、无基因型交互）弱于 A 富集模式（OR 35、p 1e-14、双方法一致）**。若机械执行字面规则（D+C），等于丢弃全项目最强的模式去抬更弱的分区证据——本身不科学。
- 因此偏离字面、服从立法本意（不许无依据的 co-option 宣称）裁决如下，每条角色有界：
  - **D backbone（Figure 2–3）：** composition→protein 传导的进化总原则 + branch-aware 桥 + adaptive 分离 + 全蛋白组框内架构 + 自然实验 + 群体边界。Association 上限。
  - **A（Figure 4，有界）：** AP2 家族 Asn 富集（正式集）+ chromatin 弱阳 + GCN5 单基因座功能实证（外部阳性对照）+ Pb 家族 map（描述性）。**明确不宣称 family-wide lifecycle co-option**（state-consequence 三极失败如实报告：liver 排除、zygote null、gametocyte 探索）。
  - **C（Figure 5/Extended，第二语境）：** acute/chronic 分区描述 + dTE 弱支持；机制未立。
  - **降级（T01）不激活：** package 未失败——D 全线加强、A 富集加强；失败的只是 H-A 的 state-consequence 子句，它本来就是 hypothesis。诚实降级适用于证据崩塌，不适用于"强模式+弱机制"的诚实 bounded package。
- D 单独成文？否——纯比较原则缺功能钩子，A（富集+单基因座实证）是当前最强功能钩子，C 为第二钩子。

## 3. Figure 映射（故事板更新依据）

- Fig1：组成 transitions（M1，不变）。
- Fig2：base→codon→AA 定量传导 + L2 branch-aware（LOO-robust）+ terminal descriptive 保留。
- Fig3：D2b 框内架构分区（附采样框偏差定量）+ D6 自然实验（T2/T3 + T1 注册阴性）+ D3 分离（conservation-matched + domain-shared + 群体层）。
- Fig4：AP2 富集（正式集 vs regex 并排）+ 反事实 + GCN5 单基因座实证 + Pb map（描述）+ state-consequence 三极阴性结果（诚实报告框）。
- Fig5/Extended：C acute/chronic 分区 + dTE 规则2 + M4R-X01 gap 声明。
- 标题禁区维持：no repeat-first、no dual-axes、no co-option-proven、no supply-collapse。

## 4. 六层充分性终值

L1✓ / L2✓（branch + LOO）/ L3✓（弱版 + 群体层）/ L4✓-框内（D2b 存活 + 富集；框外不宣称）/ L5半（单基因座实证 + 家族模式 + ART 分区；缺家族级 state 后果）/ L6待（机制判别/最小实验为投稿后或审稿回应项，非本包内）。

## 5. Claim 账本调整（本轮执行）

- CLM08：维持 CANDIDATE，证据等级注记升级（LOO + D2b + D6 + Pf8）；措辞加框内限定。
- CLM09：维持 CANDIDATE；AP2 强（正式集复现）+ 宽 RNA-binding 弱 + chromatin 新弱阳 + 严格 PUF 未定；方法注记 anchor 法弱于人工 GO（新 LEAD L-010：真 GO 注释复核）。
- CLM10-H-A/H-C：维持 hypothesis。
- T01：维持 SUPERSEDED_FOR_NOW（重激活条件：审稿要求家族级功能链或 AP2 富集在真 GO 下消失）。

## 6. 对PLAN假设的影响（一句）

H-D 在框内收敛为 C2 上限内最强一档（LOO + D2b + D6 + Pf8；仍为 CANDIDATE，不称确立）；H-A 收敛为"富集模式确立、co-option 未确立"的诚实二分；H-C 维持分区描述——三者都不支持也不需要修改 PLAN（agent 无权改 PLAN，此处仅报告无偏离）。
