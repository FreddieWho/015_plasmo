# Agent 可自动下载的数据清单与下载合同

**核验日期：** 2026-09-02  
**原则：** accession 和版本是事实源，网页 URL 只是入口；执行时必须再次核验。  
**范围：** 只收录公开、无需联系作者即可获取的数据。数据按证据门分级，不允许在项目启动时全部下载。

## 1. 存储结构

```text
data/
  raw/<source>/<dataset_id>/          # 原始下载，只读
  metadata/                           # 样本、物种、accession、ID 映射
  checksums/                          # SHA256/官方 MD5
  manual_inbox/                       # 用户人工下载的原始文件
  interim/                            # 可重建中间文件
  derived/<wp>/<analysis_id>/         # 分析产物
```

每个下载目录必须有 `DOWNLOAD_MANIFEST.json`，至少记录：
- source、accession/query、resolved URL；
- 下载时间、工具和版本；
- 文件大小与 checksum；
- 数据许可/引用；
- 是否完整、是否经过筛选；
- 对应 Gate 和使用边界。

## 2. 自动下载预算

| 预计下载量 | 默认行为 |
|---|---|
| ≤20 GB | 可自动执行，先确认剩余磁盘 |
| 20–100 GB | 仅在当前 Gate 必需且有磁盘/时间预算时执行 |
| >100 GB | 必须暂停并请求用户批准 |

**硬规则：** 优先下载 processed tables/matrices；除非重分析原始数据能区分一个关键竞争机制，否则不得下载全量 raw reads。Pf8/Pv4 等群体项目默认使用 VCF/Zarr/summary，不下载数万样本 FASTQ。

## 3. A 级：M0–M2 必需底座

### A1｜候选物种参考基因组、CDS、蛋白与注释

**主来源：** NCBI Datasets  
**入口：** https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/  
**对象：** `manifests/species_panel_seed.tsv` 中的候选物种。

**工作流：**
1. 先按 taxon 查询所有当前 assemblies；
2. 根据 assembly level、annotation、N50、完整度、参考品系与文献使用情况选择；
3. 将 assembly accession/version 写回 species manifest 并冻结；
4. 再按 accession 下载匹配的 genome/CDS/protein/GFF3/GBFF/seq-report。

示意命令（执行前以当前 `datasets --help` 核对参数）：

```bash
datasets download genome accession <GCA_or_GCF>   --include genome,cds,protein,gff3,gbff,seq-report   --filename <species>_<assembly>.zip
```

**验收：** sequence IDs 与 GFF 一致；CDS 可被 3 整除（排除标记异常者）；蛋白/CDS 数量可解释；核与细胞器分开；MD5/SHA256 保存。

### A2｜PlasmoDB/VEuPathDB 注释快照

**来源：** PlasmoDB  
**入口：** https://plasmodb.org/  
**批量目录根：** https://plasmodb.org/common/downloads/

**用途：** 寄生虫特异 gene ID、GO、产品名、定位、文献注释与 NCBI 版本交叉检查。  
**下载策略：** Agent 自动解析当前 release 目录，只下载已锁定物种的 genome/GFF/CDS/protein/annotation tables；记录 release number，不使用无版本网页结果作为永久输入。

**注意：** NCBI 与 PlasmoDB 注释有冲突时，主分析不自动择优；先做映射与差异审计，关键候选交给寄生虫专家复核。

### A3｜分类学与宿主基础元数据

**来源：** NCBI Taxonomy、权威物种论文、PlasmoDB。  
**自动对象：** taxonomy ID、标准物种名、参考品系、主要脊椎宿主、主要红细胞生态、生命周期位置。  
**边界：** 生态元数据必须带来源和置信等级，不从搜索摘要自动填充为“事实”。

### A4｜宿主血红蛋白序列

**来源：** UniProt REST 或 NCBI RefSeq/Gene。  
**对象：** 锁定宿主的 HBA/HBB（必要时含多拷贝/发育型 globin）蛋白序列。  
**用途：** 仅作为第一层营养组成近似，不等于实际氨基酸通量。  
**验收：** 物种、基因型、成熟链/前体、序列版本明确；不把人类 Hb 直接套用到所有灵长类/啮齿类宿主。

### A5｜NAR 2017 mutation-accumulation 论文公开材料

**论文：** Hamilton et al., Nucleic Acids Research 45, 1889–1901 (2017), DOI: 10.1093/nar/gkw1259  
**入口：** https://pmc.ncbi.nlm.nih.gov/articles/PMC5389722/ （若 PMCID 变化，以 DOI 检索）  
**内容：** 正文、Supplementary Tables，尤其 Table S2 中的 ENA accession；随后用 ENA Browser API 按 accession 自动下载需要的公开序列/元数据。

**优先级：**
- M0/M1：补充表和已整理突变/indel 结果；
- 仅当需要重算突变谱或重复位点时：下载 raw reads。

**注意：** 不猜测 umbrella project。先从 Supplementary Table S2 解析 accession，再下载。

## 4. B 级：Gate B 选择翻译轴后必需

### B1｜iScience 2024 mRNA/tRNA 数据

**论文：** Li et al., iScience 27, 111167 (2024), DOI: 10.1016/j.isci.2024.111167  
**GEO：** GSE226632  
**入口：** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE226632  
**FTP pattern：** `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE226nnn/GSE226632/`

**先取：** series matrix、supplementary count tables、sample metadata；只有在处理后数据不足时才用 SRA raw reads。  
**用途：** 复现 mRNA/tRNA 需求、阶段/营养相关信号和方法校准。

### B2｜iScience 2024 蛋白组

**ProteomeXchange/PRIDE：** PXD056396  
**入口：** https://www.ebi.ac.uk/pride/archive/projects/PXD056396  
**推荐方式：** 通过 PRIDE project metadata/FTP 或官方 `pridepy`/REST API 列文件；先下载 processed result/quantification 和实验设计，raw MS 仅在必要时。

**PRIDE 下载说明：** https://www.ebi.ac.uk/pride/markdownpage/pridefiledownload

### B3｜Nature Microbiology 2024 tRNA 修饰—青蒿素耐药

**论文：** Small-Saunders et al., Nature Microbiology 9, 1483–1498 (2024), DOI: 10.1038/s41564-024-01664-3  
**PMC：** https://pmc.ncbi.nlm.nih.gov/articles/PMC11153160/

**公开数据：**
- Source Data Fig. 1–5 及 Extended Data：tRNA 修饰、TMT、codon counts、RSA、IC50/IC90、heat shock 等；
- PRIDE PXD043747；
- GEO GSE151189（既往转录组）。

**入口：**
- https://www.ebi.ac.uk/pride/archive/projects/PXD043747
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151189

**下载策略：** 先自动获取 PMC source-data xlsx 和 processed proteomics/transcriptomics。论文旧版正文出现 reviewer 凭证，不得记录或使用；公开 accession 已足够。

### B4｜Malaria Cell Atlas：P. falciparum

**总入口：** https://www.malariacellatlas.org/data-sets/  
**首选基础集：** https://www.malariacellatlas.org/data_set/pf-ch10x-set1/  
**预期文件：** `pf-ch10x-set1.zip`；整合 atlas 可用 `pf.zip`。

**用途：** 血期阶段表达权重。  
**边界：** 不在第一步整合全部 4.5 万细胞；先使用处理后矩阵和官方 stage metadata，按研究问题做 pseudo-bulk。

### B5｜Malaria Cell Atlas：P. knowlesi

**入口：** https://www.malariacellatlas.org/data_set/pk-ch10x-set1/  
**预期文件：** `pk-ch10x-set1.zip` 或 atlas `pk.zip`。  
**用途：** 第二谱系血期表达与 P. falciparum 对照。

### B6｜Malaria Cell Atlas：P. berghei

**全生命周期 Smart-seq2：** https://www.malariacellatlas.org/data_set/pb-ss2-set1/  
**预期文件：** `pb-ss2-set1.zip`。  
**性发育/血期增强：** https://www.malariacellatlas.org/data_set/pb-ch10x-set2/ (`pb-ch10x-set2.zip`)。

**用途：** 传播/性发育与非灵长类谱系验证。第二数据集仅在阶段信号需要时下载。

## 5. C 级：Gate C/D 功能、必需性与耐药证据

### C1｜Science 2024 系统体外药物进化

**论文：** Luth et al., Science 386, eadk9893 (2024), DOI: 10.1126/science.adk9893  
**PMC：** https://pmc.ncbi.nlm.nih.gov/articles/PMC11809290/

**优先下载：** Supplementary Data 1–7：samples、compounds、SNV/indel、CNV、targets、gene sets、featured SNVs。  
**Raw SRA：** PRJNA1022010；论文复用的既往项目包括 PRJNA385508、PRJNA504044、PRJNA560380、PRJNA299203、PRJNA253899、PRJNA226625、PRJNA308112、PRJNA315690、PRJNA167166、SRP012591。  
**Code：** https://doi.org/10.5281/zenodo.13766513

**策略：** 默认只取 Supplementary Data 和 code；只有候选位点需要重新 calling/结构变异验证时才取部分 raw reads。论文提到可向实验室索取 evolved lines，**不纳入本项目数据清单**。

### C2｜Science 2025 P. knowlesi essential genome

**论文：** Elsworth et al., Science 387, eadq6241 (2025), DOI: 10.1126/science.adq6241  
**PMC：** https://pmc.ncbi.nlm.nih.gov/articles/PMC12104972/  
**SRA：** PRJNA1116591（Tn-seq 与 lncRNA-seq）  
**Processed browser：** https://umbibio.math.umb.edu/PkEssenDB/

**策略：** 先下载文章 Supplementary Data/processed essentiality tables 或从 PkEssenDB 导出；raw Tn-seq 仅在需重新归一化或验证候选时。

### C3｜Science 2025 跨物种 supersaturation mutagenesis

**论文：** Oberstaller et al., Science 387, eadq7347 (2025), DOI: 10.1126/science.adq7347  
**PMC：** https://pmc.ncbi.nlm.nih.gov/articles/PMC12131478/

**策略：** 自动获取 Supplementary Data 与代码；从文中/补充表解析具体 SRA accession。若没有明确 umbrella accession，禁止猜测。该数据用于判断 P. falciparum 与 P. knowlesi 必需性重编程，不作为 M0–M2 依赖。

### C4｜MalariaGEN Pf8

**入口：** https://www.malariagen.net/resource/36/  
**内容：** 33,325 个 P. falciparum 样本的开放基因组变异、样本元数据、耐药 marker/CNV 等；CC-BY 4.0。  
**论文：** DOI: 10.12688/wellcomeopenres.24031.1

**策略：** 只下载候选分析所需的 VCF/Zarr/summary/CNV/metadata；不依赖已退役或可能变化的交互式 web app，不下载全量 FASTQ。

### C5｜MalariaGEN Pv4

**入口：** https://www.malariagen.net/resource/an-open-dataset-of-iplasmodium-vivax-i-genome-variation-in-1895-worldwide-samples/  
**内容：** 1,895 个 P. vivax 样本的开放变异与元数据。  
**论文：** DOI: 10.12688/wellcomeopenres.17795.1

**用途：** 高/中 GC 人疟谱系的自然变异与选择验证。默认只取 processed variation/metadata。

## 6. D 级：按具体阶段信号触发

### D1｜孢子体翻译抑制

**论文：** Lindner et al., Nature Communications 10, 4964 (2019), DOI: 10.1038/s41467-019-12936-6  
**GEO：** GSE113582  
**PRIDE：** PXD009726、PXD009727、PXD009728、PXD009729  
**触发：** 只有当候选残差明显富集于孢子体成熟、传播或翻译抑制模块时下载。

### D2｜候选结构与功能注释

- InterPro REST/API：只对锁定候选或核心 orthogroups 查询；
- AlphaFold DB/RCSB PDB：只下载候选蛋白结构；
- 禁止在 M0–M2 批量下载所有物种全部预测结构或整个 Pfam 数据库，除非有独立方法需求和预算。

## 7. 下载工具合同

### NCBI GEO/SRA
- GEO processed：FTP/HTTPS 或 GEOquery/官方 MINiML；
- SRA：先用 RunInfo 生成样本清单，`prefetch` + `fasterq-dump` 仅针对批准的 runs；
- 每个 SRR 记录 BioProject、BioSample、library strategy、layout 与 sample label。

### ENA
- 由 accession 驱动使用 ENA Browser API；
- 不通过物种关键词盲目抓取；
- 下载前确认 study/sample/run 层级。

### PRIDE
- 以 PXD accession 获取 project metadata 与 file list；
- 优先 processed/result/experimental design；
- 文件可通过 PRIDE 官方 FTP/HTTPS/streaming API 或官方 `pridepy` 获取；
- raw MS 下载前需说明为何 processed 文件不足。

### PMC supplementary data
- 优先 PMC Associated Data/Source Data 的公开附件；
- 保存原始文件名和文章 DOI；
- 自动链接失败时转入 `06_DATA_MANUAL_ACQUISITION.md`，不向作者索取。

## 8. 自动下载后验收清单

- [ ] accession、版本、来源、许可已记录；
- [ ] 文件大小非零，checksum 已生成/比对；
- [ ] 样本数量与论文/页面一致或差异已解释；
- [ ] processed/raw 状态明确；
- [ ] gene ID/assembly 版本可映射；
- [ ] 下载范围与当前 Gate 对齐；
- [ ] 没有混入作者请求、私有 reviewer 或受控数据；
- [ ] data registry 状态已更新为 `DOWNLOADED_VERIFIED` 或 `DOWNLOADED_UNRESOLVED`。
