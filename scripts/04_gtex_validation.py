import pandas as pd
import numpy as np
from pathlib import Path

BASE = Path.home() / "project/oadrm"
GTEX_FILE = BASE / "data/raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
SIG_DIR = BASE / "data/processed/organ_signatures"
PROC = BASE / "data/processed"
PROC.mkdir(parents=True, exist_ok=True)

# ==================== 1. 读取 GTEx v8 ====================
print("=== 读取 GTEx v8 基因中位数表达 ===")
print(f"文件: {GTEX_FILE}")

# GCT 格式：前两行是头信息，第三行是列名
# 第1行: #1.2
# 第2行: n_rows n_cols
# 第3行: Name Description 组织1 组织2 ...

with open(GTEX_FILE, 'r') as f:
    line1 = f.readline().strip()
    line2 = f.readline().strip()
    n_rows, n_cols = map(int, line2.split())
    print(f"头信息: {line1}, {n_rows} genes x {n_cols} tissues")

# 读取数据
gtex = pd.read_csv(GTEX_FILE, sep='\t', skiprows=2)
print(f"GTEx shape: {gtex.shape}")
print(f"Columns: {gtex.columns.tolist()[:5]} ... {gtex.columns.tolist()[-3:]}")

# 提取基因名和组织列
# 第一列是 gene_id (如 ENSG00000223972.5)，第二列是 Description
gene_col = gtex.columns[0]  # Name
desc_col = gtex.columns[1]  # Description
tissue_cols = gtex.columns[2:].tolist()

print(f"\n组织数: {len(tissue_cols)}")
print(f"组织列表: {tissue_cols}")

# 提取 Gene Symbol（从 Description 列）
gtex['gene_symbol'] = gtex[desc_col].str.split('.').str[0]
print(f"\nGene Symbol 示例: {gtex['gene_symbol'].head().tolist()}")

# ==================== 2. 器官-组织映射 ====================
# GTEx 组织名称 → 我们的器官名称
ORGAN_TISSUE_MAP = {
    'Brain': ['Brain - Cerebellum', 'Brain - Cerebellar Hemisphere', 'Brain - Cortex',
              'Brain - Frontal Cortex (BA9)', 'Brain - Anterior cingulate cortex (BA24)',
              'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Amygdala',
              'Brain - Caudate (basal ganglia)', 'Brain - Nucleus accumbens (basal ganglia)',
              'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
              'Brain - Substantia nigra'],
    'Blood': ['Whole Blood', 'Cells - EBV-transformed lymphocytes'],
    'Intestine': ['Colon - Sigmoid', 'Colon - Transverse', 'Small Intestine - Terminal Ileum'],
    'Kidney': ['Kidney - Cortex', 'Kidney - Medulla'],
    'Liver': ['Liver'],
    'Lung': ['Lung'],
    'Skin': ['Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)'],
    'Stomach': ['Stomach'],
    'Heart': ['Heart - Atrial Appendage', 'Heart - Left Ventricle'],
    'Muscle': ['Muscle - Skeletal'],
    'Pancreas': ['Pancreas'],
}

# ==================== 3. 验证每个器官签名 ====================
print("\n=== GTEx 组织表达验证 ===")

# 读取汇总表
summary = pd.read_csv(SIG_DIR / "signature_summary_v2.csv")
results = []

for _, row in summary.iterrows():
    organ = row['organ']
    if organ == 'Conventional':
        continue  # 全身签名不做GTEx验证

    # 读取签名
    sig_file = SIG_DIR / row['file']
    sig = pd.read_csv(sig_file)
    genes = sig['gene_symbol'].tolist()
    print(f"\n--- {organ} ({len(genes)} genes) ---")

    # 找到GTEx对应组织
    tissues = ORGAN_TISSUE_MAP.get(organ, [])
    if not tissues:
        print(f"  ⚠ 无GTEx组织映射")
        continue

    # 检查这些基因在GTEx中是否存在
    matched = gtex[gtex['gene_symbol'].isin(genes)]
    print(f"  GTEx匹配基因: {len(matched)}/{len(genes)}")

    if len(matched) == 0:
        print(f"  ⚠ 无基因匹配")
        continue

    # 计算各组织的中位数TPM
    tissue_data = matched[tissues + ['gene_symbol']].set_index('gene_symbol')

    # 计算该器官组织 vs 所有组织的表达排名
    all_tissue_median = gtex[tissue_cols].median(axis=1)

    validation = []
    for gene in genes:
        if gene not in tissue_data.index:
            continue

        # 该基因在靶器官的表达（取多个组织的最大值）
        gene_organ_expr = tissue_data.loc[gene].max()

        # 该基因在所有组织中的中位数表达
        gene_all_median = all_tissue_median[gtex['gene_symbol'] == gene].values
        if len(gene_all_median) == 0:
            continue
        gene_all_median = gene_all_median[0]

        # 计算百分位排名（在所有组织中）
        gene_all_values = gtex[gtex['gene_symbol'] == gene][tissue_cols].values.flatten()
        if len(gene_all_values) == 0:
            continue
        percentile = (gene_organ_expr >= gene_all_values).mean() * 100

        validation.append({
            'gene_symbol': gene,
            'organ_expr': gene_organ_expr,
            'all_median': gene_all_median,
            'percentile': percentile,
            'pass': percentile >= 75  # 前25%阈值
        })

    val_df = pd.DataFrame(validation)
    if len(val_df) > 0:
        pass_rate = val_df['pass'].mean() * 100
        median_pct = val_df['percentile'].median()
        print(f"  验证基因: {len(val_df)}")
        print(f"  中位数百分位: {median_pct:.1f}%")
        print(f"  通过率(≥75%): {pass_rate:.1f}%")

        results.append({
            'organ': organ,
            'n_genes': len(genes),
            'n_matched': len(matched),
            'n_validated': len(val_df),
            'median_percentile': median_pct,
            'pass_rate': pass_rate
        })

        # 保存详细验证结果
        val_df.to_csv(PROC / f"gtex_validation_{organ.lower()}.csv", index=False)
    else:
        print(f"  ⚠ 无验证结果")

# ==================== 4. 保存汇总 ====================
print("\n=== GTEx 验证汇总 ===")
if results:
    result_df = pd.DataFrame(results)
    result_df.to_csv(PROC / "gtex_validation_summary.csv", index=False)
    print(result_df.to_string(index=False))
else:
    print("无验证结果")

print("\n=== 完成 ===")
