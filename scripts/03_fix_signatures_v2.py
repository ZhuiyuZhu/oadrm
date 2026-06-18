import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

BASE = Path.home() / "project/oadrm"
RAW = BASE / "data/raw/nc_supp"
PROC = BASE / "data/processed/organ_signatures"
PROC.mkdir(parents=True, exist_ok=True)

# ==================== 1. 重新读取 Data 3 ====================
print("=== 重新读取 Supplementary Data 3 ===")
FILE = RAW / "41467_2025_67223_MOESM2_ESM.xlsx"

df = pd.read_excel(FILE, sheet_name='Supplementary Data 3',
                   header=None, skiprows=2,
                   names=['organ', 'gene_symbol', 'weight'])

df['organ'] = df['organ'].ffill()
df['weight'] = pd.to_numeric(df['weight'], errors='coerce')
df = df.dropna()
df = df[df['gene_symbol'] != 'sex']
df = df[df['gene_symbol'].str.len() >= 2]

print(f"清洗后总记录: {len(df)}")
print(f"\n器官分布:\n{df['organ'].value_counts()}")

# ==================== 2. 器官映射 ====================
ORGAN_MAP = {
    'Brain': 'Brain', 'Intestine': 'Intestine', 'Kidney': 'Kidney',
    'Liver': 'Liver', 'Lung': 'Lung', 'Skin': 'Skin', 'Stomach': 'Stomach',
    'Immune': 'Blood', 'Conventional': 'Conventional',
    'Heart': 'Heart', 'Muscle': 'Muscle', 'Pancreas': 'Pancreas',
}

df['target_organ'] = df['organ'].map(ORGAN_MAP)
df = df.dropna(subset=['target_organ'])

LINCS_AVAILABLE = ['Brain', 'Intestine', 'Kidney', 'Liver', 'Lung',
                   'Skin', 'Stomach', 'Blood']
available = [o for o in sorted(df['target_organ'].unique())
             if o in LINCS_AVAILABLE or o == 'Conventional']
print(f"\n✓ 核心可用器官 ({len(available)} 个): {available}")

# ==================== 3. 提取签名并保存 ====================
print("\n=== 提取并保存签名 ===")
summary = []

for organ in available:
    organ_df = df[df['target_organ'] == organ].copy()
    organ_df['abs_weight'] = organ_df['weight'].abs()

    n_top = min(100, len(organ_df))
    top_df = organ_df.nlargest(n_top, 'abs_weight').copy()

    max_w = top_df['weight'].abs().max()
    top_df['w_norm'] = top_df['weight'] / max_w if max_w > 0 else 0
    top_df['direction'] = np.sign(top_df['weight'])

    # 修复：先算绝对值，再 groupby 取 idxmax
    top_df['abs_w'] = top_df['weight'].abs()
    idx = top_df.groupby('gene_symbol')['abs_w'].idxmax()
    top_df = top_df.loc[idx].copy()

    out_file = PROC / f"signature_{organ.lower()}_top{len(top_df)}.csv"
    top_df[['gene_symbol', 'weight', 'w_norm', 'direction']].to_csv(out_file, index=False)

    n_up = sum(top_df['direction'] > 0)
    n_down = sum(top_df['direction'] < 0)
    summary.append({
        'organ': organ, 'n_genes': len(top_df),
        'n_up': n_up, 'n_down': n_down, 'file': out_file.name
    })
    print(f"  ✓ {organ:15s}: {len(top_df):3d} genes (up={n_up}, down={n_down}) → {out_file.name}")

summary_df = pd.DataFrame(summary)
summary_df.to_csv(PROC / "signature_summary_v2.csv", index=False)
print(f"\n✓ 汇总表: {PROC / 'signature_summary_v2.csv'}")
print(summary_df.to_string(index=False))

# ==================== 4. Bootstrap稳健性 ====================
print("\n=== Bootstrap稳健性评估 ===")
np.random.seed(42)

for organ in available:
    if organ == 'Conventional':
        continue
    organ_df = df[df['target_organ'] == organ].copy()
    n_boot = 100
    selected_sets = []

    for _ in range(n_boot):
        boot = organ_df.sample(frac=0.8, replace=True)
        boot['abs_w'] = boot['weight'].abs()
        n = min(100, len(boot))
        top = set(boot.nlargest(n, 'abs_w')['gene_symbol'])
        selected_sets.append(top)

    gene_freq = Counter([g for s in selected_sets for g in s])
    stable = sum(1 for f in gene_freq.values() if f >= 0.8 * n_boot)
    print(f"  {organ:15s}: {stable}/{n_boot} 基因稳定 ({stable}%)")

print("\n=== 完成 ===")
