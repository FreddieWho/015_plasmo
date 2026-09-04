# 可能需要用户人工操作的数据清单

**核验日期：** 2026-09-02  
**核心结论：** 当前方案没有任何“必须人工申请或联系作者”才能启动的数据。人工操作只作为公开数据自动获取失败时的后备路径。

## 1. 用户操作统一规则

1. 只从本文件列出的官方网页/仓库下载；
2. 保留原始文件名，不重命名、不解压覆盖；
3. 放入 `data/manual_inbox/<item_id>/`；
4. 记录下载日期、页面 URL、文件名和页面显示的版本；
5. 生成 SHA256：

```bash
sha256sum <file> > <file>.sha256
```

6. 通知 Agent 后，由 Agent 完成内容核验和登记；用户不需要自行整理表格。

## 2. 人工后备队列

### M01｜NAR 2017 Supplementary Table S2

**何时触发：** Agent 无法通过 DOI/PMC/publisher 自动取得 supplement，因而无法解析 ENA accession。  
**官方入口：** DOI 10.1093/nar/gkw1259；NAR/PMC 文章页面。  
**需要下载：** 完整 Supplementary Data，尤其包含 accession 的 Table S2。  
**目标目录：** `data/manual_inbox/M01_NAR2017_SUPPLEMENT/`  
**说明：** raw 序列仍由 Agent 根据公开 ENA accession 自动下载；不需要联系作者。

### M02｜当前 PlasmoDB release 的物种批量文件

**何时触发：** 官方版本目录可浏览，但自动脚本无法解析当前 release 或被动态页面/反爬机制阻断。  
**官方入口：** https://plasmodb.org/ → Downloads；或 https://plasmodb.org/common/downloads/  
**需要下载：** 仅已锁定物种的 genome、GFF、CDS、protein 和 annotation table。  
**必须同时记录：** release number、物种/品系名、每个文件网页路径。  
**目标目录：** `data/manual_inbox/M02_PLASMODB_RELEASE_<N>/`

### M03｜Nature Microbiology 2024 Source Data

**何时触发：** PMC 自动附件下载失败。  
**官方入口：** https://pmc.ncbi.nlm.nih.gov/articles/PMC11153160/  
**需要下载：** Source Data Fig. 1–5、相关 Extended Data xlsx、Supplementary Information；不需要下载 `.ai` 原始图文件，除非图像审计明确要求。  
**目标目录：** `data/manual_inbox/M03_NATMICO2024_SOURCE_DATA/`

### M04｜Science 2024 Luth Supplementary Data

**何时触发：** PMC 自动附件下载失败。  
**官方入口：** https://pmc.ncbi.nlm.nih.gov/articles/PMC11809290/  
**需要下载：** Supplementary Data 1–7 与 Supplemental Material PDF。  
**目标目录：** `data/manual_inbox/M04_LUTH2024_SUPPLEMENT/`

### M05｜Malaria Cell Atlas 处理后数据

**何时触发：** 下载链接由动态页面生成，Agent 无法取得 zip。  
**官方入口：** https://www.malariacellatlas.org/data-sets/  
**按当前计划只需：**
- `pf-ch10x-set1.zip`；
- `pk-ch10x-set1.zip`；
- `pb-ss2-set1.zip`；
- `pb-ch10x-set2.zip` 仅条件触发。

**目标目录：** `data/manual_inbox/M05_MALARIA_CELL_ATLAS/`

### M06｜MalariaGEN 开放数据包/条款点击

**何时触发：** 数据仍是开放数据，但下载入口迁移、临时签名 URL 或网页条款需要浏览器确认。  
**官方入口：** Pf8 https://www.malariagen.net/resource/36/；Pv4 官方 resource 页面。  
**需要下载：** processed VCF/Zarr/metadata/CNV 或页面提供的 release bundle，不下载全量 raw FASTQ。  
**目标目录：** `data/manual_inbox/M06_MALARIAGEN/`  
**触发前提：** 当前 Gate 明确需要群体证据，且 Agent 已给出文件规模。

## 3. 明确不纳入的资料

以下数据即使在论文中被提及，也不得进入下载队列：

- 需要给作者发邮件、填写“合理请求”或等待实验室批准的数据；
- Luth et al. 论文中需向 originating lab 索取的 evolved parasite lines；
- iScience 论文中仅称“additional information available from lead contact”的未公开材料；
- 未公开、manuscript in preparation 且没有公开下载文件的 Malaria Cell Atlas 数据；
- 需要作者单独授权的私有 Zenodo/Drive/Dropbox 文件；
- EGA/dbGaP 等受控人类宿主组学和患者级数据，除非未来形成独立且必要的获批子课题；
- 通过 reviewer username/password 才能访问的数据；若 accession 已公开，只使用公开入口；
- 任何来历不明的二次网盘镜像。

## 4. 人工数据验收

Agent 收到人工文件后必须：
- 验证文件与论文/项目相符；
- 检查是否公开及许可；
- 核对 checksum、文件数、样本数；
- 写入 `data_registry.tsv`；
- 将原文件移入只读 raw 目录的副本，保留 manual_inbox 原件；
- 若仍缺失，不得请求作者，必须选择替代公开证据、调整分析或记录为不可得。
