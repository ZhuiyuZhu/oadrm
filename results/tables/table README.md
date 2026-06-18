# results/tables/ — 结果表格与矩阵数据

## 中文

本目录存放 OADRM 筛选和统计产生的**结构化表格**（CSV 格式）。

| 文件名 | 内容说明 | 大小/备注 |
|-------|---------|----------|
| `organ_drug_top20_v3_cleaned.csv` | 每个器官 Top 20 候选药物（过滤后） | 小文件，可直接下载 |
| `organ_drug_full_matrix_v3.csv` | 完整 OADRM 矩阵：8 器官 × 21,304 化合物 | **~500MB–1GB**，通常不上传 GitHub |
| `homologous_vs_heterologous_stats.csv` | 同源 vs. 异源 mCS 的统计检验结果（Wilcoxon p 值、Cohen's d） | 小文件 |

### ⚠️ 大文件处理
`organ_drug_full_matrix_v3.csv` 为完整评分矩阵，体积过大，**默认不纳入 Git 版本控制**。请通过以下途径获取：
- **Zenodo**: [10.5281/zenodo.xxxxxxx](https://doi.org/10.5281/zenodo.xxxxxxx)
- **Figshare**: [xxxxxxx](https://figshare.com/xxxxxxx)
- 或运行 `scripts/06_connectivity_score_lowmem.py` + `scripts/07_drug_prioritization.py` 自行生成

---

## English

This directory stores **structured tables** (CSV) produced by OADRM screening and statistical analysis.

| Filename | Description | Size / Note |
|---------|-------------|-------------|
| `organ_drug_top20_v3_cleaned.csv` | Top 20 candidates per organ (after filtering) | Small file, downloadable directly |
| `organ_drug_full_matrix_v3.csv` | Full OADRM matrix: 8 organs × 21,304 compounds | **~500MB–1GB**, usually not uploaded to GitHub |
| `homologous_vs_heterologous_stats.csv` | Statistical test results for homologous vs. heterologous mCS (Wilcoxon p, Cohen's d) | Small file |

### ⚠️ Large File Handling
`organ_drug_full_matrix_v3.csv` is the complete scoring matrix and is too large for Git. Obtain it via:
- **Zenodo**: [10.5281/zenodo.xxxxxxx](https://doi.org/10.5281/zenodo.xxxxxxx)
- **Figshare**: [xxxxxxx](https://figshare.com/xxxxxxx)
- Or generate locally by running `scripts/06_connectivity_score_lowmem.py` + `scripts/07_drug_prioritization.py`
