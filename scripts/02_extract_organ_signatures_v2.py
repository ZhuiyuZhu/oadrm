import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

BASE = Path.home() / "project/oadrm"
RAW = BASE / "data/raw/nc_supp"
PROC = BASE / "data/processed/organ_signatures"
PROC.mkdir(parents=True, exist_ok=True)

# ==================== 1. 读取 Data 3 ====================
print("=== 读取 Supplementary Data 3 ===")
FILE = RAW / "41467_2025_67223_MOESM2_ESM.xlsx"

df = pd.read_excel(FILE, sheet_name='Supplementary Data 3',
                   header=None, skiprows=2,
                   names=['organ', 'protein_id', 'weight'])

df['organ'] = df['organ'].ffill()
df['weight'] = pd.to_numeric(df['weight'], errors='coerce')
df = df.dropna()

print(f"总记录数: {len(df)}")
print(f"\n器官分布:\n{df['organ'].value_counts()}")

# ==================== 2. 器官映射策略 ====================
# Data 3 中的器官 → 目标器官名称
# 策略：保留与LINCS直接匹配的，可选映射近似的

ORGAN_MAP = {
    # 直接匹配（LINCS有对应细胞系）
    'Brain': 'Brain',
    'Intestine': 'Intestine',
    'Kidney': 'Kidney',
    'Liver': 'Liver',
    'Lung': 'Lung',
    'Skin': 'Skin',
    'Stomach': 'Stomach',
    # 近似映射（可选，如不需要可注释掉）
    'Immune': 'Blood',  # 免疫系统 → 血液细胞系
    'Heart': 'Heart',  # LINCS中无心肌细胞系，但保留供后续补充
    'Muscle': 'Muscle',  # LINCS中无骨骼肌细胞系，保留
    'Pancreas': 'Pancreas',  # LINCS中无胰腺细胞系，保留
    # Conventional 作为全身/广谱参考（单独处理）
    'Conventional': 'Conventional',
    # 以下器官LINCS无直接对应，暂排除
    # 'Adipose': None,
    # 'Artery': None,
}

# 应用映射
df['target_organ'] = df['organ'].map(ORGAN_MAP)
df_mapped = df.dropna(subset=['target_organ']).copy()

print(f"\n=== 映射后可用器官 ===")
available_organs = sorted(df_mapped['target_organ'].unique())
print(f"可用器官 ({len(available_organs)} 个): {available_organs}")

for organ in available_organs:
    count = len(df_mapped[df_mapped['target_organ'] == organ])
    print(f"  {organ:15s}: {count:4d} 个蛋白")

# ==================== 3. 提取各器官Top 100签名 ====================
print("\n=== 提取器官签名 ===")
signatures = {}

for organ in available_organs:
    organ_df = df_mapped[df_mapped['target_organ'] == organ].copy()

    # 按绝对权重取Top 100
    organ_df['abs_weight'] = organ_df['weight'].abs()
    top_df = organ_df.nlargest(min(100, len(organ_df)), 'abs_weight').copy()

    # 标准化权重
    max_abs_w = top_df['weight'].abs().max()
    top_df['w_norm'] = top_df['weight'] / max_abs_w if max_abs_w > 0 else 0

    # 方向：+1=随龄上调，-1=下调
    top_df['direction'] = np.sign(top_df['weight'])

    signatures[organ] = top_df[['protein_id', 'weight', 'w_norm', 'direction']]

    n_up = sum(top_df['direction'] > 0)
    n_down = sum(top_df['direction'] < 0)  # 修复：<< 改为 <
    print(f"  {organ:15s}: {len(top_df):3d} proteins, up={n_up}, down={n_down}")

# ==================== 4. 蛋白→Gene Symbol映射 ====================
print("\n=== 蛋白→Gene Symbol映射 ===")

try:
    import mygene

    mg = mygene.MyGeneInfo()

    all_proteins = set()
    for organ, sig_df in signatures.items():
        all_proteins.update(sig_df['protein_id'].astype(str).tolist())

    all_proteins = list(all_proteins)
    print(f"待映射蛋白数: {len(all_proteins)}")

    # 判断ID格式
    sample = all_proteins[0]
    print(f"蛋白ID示例: {sample} (长度{len(sample)})")

    # 批量查询
    batch_size = 1000
    mapping_results = []

    for i in range(0, len(all_proteins), batch_size):
        batch = all_proteins[i:i + batch_size]
        result = mg.querymany(batch, scopes='uniprot',
                              fields='symbol,name,ensembl.gene',
                              species='human', verbose=False)
        mapping_results.extend(result)
        if (i // batch_size + 1) % 5 == 0 or i == 0:
            print(f"  进度: {min(i + batch_size, len(all_proteins))}/{len(all_proteins)}")

    # 构建映射表
    map_df = pd.DataFrame(mapping_results)
    if 'query' in map_df.columns and 'symbol' in map_df.columns:
        map_df = map_df[['query', 'symbol', 'ensembl']].drop_duplicates('query')
        map_df.columns = ['protein_id', 'gene_symbol', 'ensembl_id']

        mapped_signatures = {}
        for organ, sig_df in signatures.items():
            merged = sig_df.merge(map_df, on='protein_id', how='left')
            mapped_signatures[organ] = merged
            mapped = merged['gene_symbol'].notna().sum()
            print(f"  {organ:15s}: {mapped}/{len(merged)} 映射成功")
    else:
        print("⚠ mygene返回格式异常")
        mapped_signatures = signatures

except ImportError:
    print("⚠ mygene未安装，运行: pip install mygene")
    mapped_signatures = signatures
except Exception as e:
    print(f"⚠ 映射失败: {e}")
    mapped_signatures = signatures

# ==================== 5. 保存签名文件 ====================
print("\n=== 保存签名 ===")
summary = []

for organ, sig_df in mapped_signatures.items():
    if 'gene_symbol' in sig_df.columns:
        out_df = sig_df.dropna(subset=['gene_symbol']).copy()
        if len(out_df) > 0:
            out_df = out_df.loc[out_df.groupby('gene_symbol')['weight'].abs().idxmax()]
        out_file = PROC / f"signature_{organ.lower().replace(' ', '_')}_top100.csv"
        out_df[['gene_symbol', 'weight', 'w_norm', 'direction', 'protein_id']].to_csv(
            out_file, index=False
        )
        n_genes = len(out_df)
    else:
        out_file = PROC / f"signature_{organ.lower().replace(' ', '_')}_top100_protein.csv"
        sig_df.to_csv(out_file, index=False)
        n_genes = len(sig_df)

    summary.append({
        'organ': organ,
        'n_proteins': len(sig_df),
        'n_genes_mapped': n_genes,
        'n_up': sum(sig_df['direction'] > 0),
        'n_down': sum(sig_df['direction'] < 0)
    })
    print(f"  ✓ {organ:15s}: {n_genes} → {out_file.name}")

# 保存汇总
summary_df = pd.DataFrame(summary)
summary_df.to_csv(PROC / "signature_summary.csv", index=False)
print(f"\n✓ 汇总表: {PROC / 'signature_summary.csv'}")
print(summary_df.to_string(index=False))

# ==================== 6. Bootstrap稳健性（示例） ====================
print("\n=== Bootstrap稳健性评估（前3个器官示例） ===")
np.random.seed(42)

for organ in list(signatures.keys())[:3]:
    organ_df = df_mapped[df_mapped['target_organ'] == organ].copy()
    n_boot = 100
    selected_sets = []

    for _ in range(n_boot):
        boot = organ_df.sample(frac=0.8, replace=True)
        boot['abs_w'] = boot['weight'].abs()
        top_genes = set(boot.nlargest(min(100, len(boot)), 'abs_w')['protein_id'].astype(str))
        selected_sets.append(top_genes)

    gene_freq = Counter([g for s in selected_sets for g in s])
    stable = sum(1 for f in gene_freq.values() if f >= 0.8 * n_boot)
    print(f"  {organ:15s}: {stable}/100 基因稳定 ({stable}%)")

print("\n=== 完成 ===")
