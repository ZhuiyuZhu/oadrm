# data/ — 输入数据与中间处理结果

## 中文

本目录存放 OADRM 流水线所需的**输入数据**和**中间处理结果**。

| 子目录/文件 | 内容说明 | 备注 |
|-----------|---------|------|
| `processed/organ_signatures/` | 8 个器官的衰老转录组特征（CSV 格式） | 从 UKB-PPP 蛋白时钟提取，经 UniProt 映射和 GTEx v8 组织富集过滤后生成 |
| `processed/lincs_meta/` | LINCS L1000 细胞系-器官映射表 | 83 个细胞系手动标注的组织来源（如 HepG2→liver, A549→lung） |
| `processed/gtex_validation_summary.csv` | GTEx v8 组织富集验证汇总 | 各器官特征基因在目标组织中的 TPM 排名百分位及通过率 |

### ⚠️ 重要声明
- **本仓库不存放 UK Biobank 原始个体级数据**（受 UKB 数据使用协议保护）。
- **本仓库不存放 LINCS L1000 原始 Level-5 矩阵**（文件过大，请从 GEO [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) 和 [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138) 自行下载）。
- 仅上传**小体积的衍生结果文件**（特征列表、映射表、验证统计），以便复现和审计。

---

## English

This directory stores **input data** and **intermediate processed results** for the OADRM pipeline.

| Subdirectory/File | Description | Note |
|------------------|-------------|------|
| `processed/organ_signatures/` | Transcriptional aging signatures for 8 organs (CSV) | Derived from UKB-PPP protein clocks, mapped via UniProt, and filtered by GTEx v8 tissue enrichment |
| `processed/lincs_meta/` | LINCS L1000 cell line–organ mapping table | Manual annotation of tissue origin for 83 cell lines (e.g., HepG2→liver, A549→lung) |
| `processed/gtex_validation_summary.csv` | GTEx v8 tissue enrichment validation summary | Percentile rank and pass rate of signature genes in target tissues |

### ⚠️ Important Notes
- **Raw UK Biobank individual-level data are NOT stored here** (protected by UKB data agreement).
- **Raw LINCS L1000 Level-5 matrices are NOT stored here** (too large; download from GEO [GSE92742](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742) & [GSE70138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138)).
- Only **small derivative files** (signature lists, mapping tables, validation statistics) are uploaded for reproducibility and audit.
