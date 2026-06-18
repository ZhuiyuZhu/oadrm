import pandas as pd
import numpy as np
from pathlib import Path
import glob

BASE = Path.home() / "project/oadrm"
GTEX_FILE = BASE / "data/raw/gtex/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
SIG_DIR = BASE / "data/processed/organ_signatures"
PROC = BASE / "data/processed"
PROC.mkdir(parents=True, exist_ok=True)

# ==================== 1. LINCS 基因匹配 ====================
print("=== LINCS 基因匹配 ===")
g1 = pd.read_csv(BASE / "data/raw/lincs/GSE92742_Broad_LINCS_gene_info.txt", sep='\t')
g2 = pd.read_csv(BASE / "data/raw/lincs/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", sep='\t')
lincs_genes = pd.concat([g1, g2]).drop_duplicates(subset='pr_gene_symbol')
lincs_gene_set = set(lincs_genes['pr_gene_symbol'].dropna().astype(str).str.upper())

print(f"LINCS总基因数: {len(lincs_gene_set)}")

sig_files = sorted(glob.glob(str(SIG_DIR / "signature_*_top*.csv")))
sig_files = [f for f in sig_files if '_protein' not in f and 'summary' not in f]
lincs_results = []

for f in sig_files:
    fname = Path(f).name
    sig = pd.read_csv(f)
    organ = fname.replace('signature_', '').replace('.csv', '')
    our_genes = set(sig['gene_symbol'].dropna().astype(str).str.upper())
    matched = our_genes & lincs_gene_set
    pct = len(matched) / len(our_genes) * 100 if our_genes else 0

    lincs_results.append({
        'organ': organ, 'n_sig': len(our_genes),
        'n_matched': len(matched), 'pct': pct
    })
    print(f"  {organ:25s}: {len(our_genes):3d}个, 匹配{len(matched):3d}个 ({pct:5.1f}%)")

lincs_df = pd.DataFrame(lincs_results)
lincs_df.to_csv(PROC / "lincs_gene_match_summary.csv", index=False)
print(f"\n✓ LINCS匹配汇总: {PROC / 'lincs_gene_match_summary.csv'}")

# ==================== 2. GTEx 组织表达验证 ====================
print("\n=== GTEx 组织表达验证 ===")

# 正确读取：跳过前2行（#1.2 和 56200 54），第3行是表头
gtex = pd.read_csv(GTEX_FILE, sep='\t', skiprows=2)
print(f"GTEx shape: {gtex.shape}")
print(f"Columns前5: {gtex.columns.tolist()[:5]}")

desc_col = gtex.columns[1]  # Description
tissue_cols = gtex.columns[2:].tolist()
gtex['gene_symbol'] = gtex[desc_col].str.split('.').str[0].str.upper()


def clean_organ(name):
    return name.replace('_top100', '').replace('_top79', '').replace('_top34', '').replace('_top42', '').replace(
        '_top70', '').replace('_top38', '').replace('_top96', '')


ORGAN_TISSUE = {
    'blood': ['Whole Blood', 'Cells - EBV-transformed lymphocytes'],
    'brain': ['Brain - Cerebellum', 'Brain - Cerebellar Hemisphere', 'Brain - Cortex',
              'Brain - Frontal Cortex (BA9)', 'Brain - Anterior cingulate cortex (BA24)',
              'Brain - Hippocampus', 'Brain - Hypothalamus', 'Brain - Amygdala',
              'Brain - Caudate (basal ganglia)', 'Brain - Nucleus accumbens (basal ganglia)',
              'Brain - Putamen (basal ganglia)', 'Brain - Spinal cord (cervical c-1)',
              'Brain - Substantia nigra'],
    'intestine': ['Colon - Sigmoid', 'Colon - Transverse', 'Small Intestine - Terminal Ileum'],
    'kidney': ['Kidney - Cortex', 'Kidney - Medulla'],
    'liver': ['Liver'],
    'lung': ['Lung'],
    'skin': ['Skin - Not Sun Exposed (Suprapubic)', 'Skin - Sun Exposed (Lower leg)'],
    'stomach': ['Stomach'],
}

gtex_results = []

for _, row in lincs_df.iterrows():
    organ = row['organ']
    organ_clean = clean_organ(organ)

    if organ_clean == 'conventional':
        continue

    f = SIG_DIR / f"signature_{organ}.csv"
    sig = pd.read_csv(f)
    genes = sig['gene_symbol'].dropna().astype(str).str.upper().tolist()

    tissues = ORGAN_TISSUE.get(organ_clean, [])
    if not tissues:
        print(f"\n--- {organ_clean}: 无GTEx组织映射 ---")
        continue

    print(f"\n--- {organ_clean} ({len(genes)} genes) ---")

    matched = gtex[gtex['gene_symbol'].isin(genes)]
    print(f"  GTEx匹配: {len(matched)}/{len(genes)}")

    if len(matched) == 0:
        continue

    validations = []
    for gene in genes:
        gene_rows = matched[matched['gene_symbol'] == gene]
        if len(gene_rows) == 0:
            continue

        organ_expr_vals = gene_rows[tissues].values.flatten()
        organ_max = np.nanmax(organ_expr_vals) if len(organ_expr_vals) > 0 else 0

        all_expr = gene_rows[tissue_cols].values.flatten()
        all_median = np.nanmedian(all_expr)

        percentile = np.mean(organ_max >= all_expr) * 100 if len(all_expr) > 0 else 0

        validations.append({
            'gene_symbol': gene,
            'organ_max_tpm': organ_max,
            'all_median_tpm': all_median,
            'percentile': percentile,
            'pass': percentile >= 75
        })

    val_df = pd.DataFrame(validations)
    if len(val_df) > 0:
        pass_rate = val_df['pass'].mean() * 100
        median_pct = val_df['percentile'].median()
        print(f"  验证基因: {len(val_df)}")
        print(f"  中位数百分位: {median_pct:.1f}%")
        print(f"  通过率(≥75%): {pass_rate:.1f}%")

        gtex_results.append({
            'organ': organ_clean, 'n_sig': len(genes), 'n_matched': len(matched),
            'n_validated': len(val_df), 'median_pct': median_pct, 'pass_rate': pass_rate
        })

        val_df.to_csv(PROC / f"gtex_validation_{organ_clean}.csv", index=False)

if gtex_results:
    gtex_summary = pd.DataFrame(gtex_results)
    gtex_summary.to_csv(PROC / "gtex_validation_summary.csv", index=False)
    print(f"\n=== GTEx 验证汇总 ===")
    print(gtex_summary.to_string(index=False))
    print(f"\n✓ GTEx汇总: {PROC / 'gtex_validation_summary.csv'}")

print("\n=== 完成 ===")
