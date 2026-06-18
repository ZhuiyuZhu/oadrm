import pandas as pd
import numpy as np
from pathlib import Path
import glob
from scipy import stats

BASE = Path.home() / "project/oadrm"
SCORE_DIR = BASE / "data/processed/drug_scores"
META_DIR = BASE / "data/processed/lincs_meta"
RESULTS = BASE / "results/tables"
RESULTS.mkdir(parents=True, exist_ok=True)

# 1. 读取所有 mCS 结果
print("=== 读取 mCS 数据 ===")
mcs_files = sorted(glob.glob(str(SCORE_DIR / "*_mcs.csv")))
all_scores = []
for f in mcs_files:
    fname = Path(f).stem
    parts = fname.split('_')
    dataset, organ = parts[0], parts[1]
    df = pd.read_csv(f)
    df['dataset'] = dataset
    df['target_organ'] = organ.lower()
    all_scores.append(df)

scores = pd.concat(all_scores, ignore_index=True)
print(f"合并: {len(scores)} rows")

# 2. 读取细胞系-器官映射
cell_organ = pd.read_csv(META_DIR / "cellline_to_organ.csv")
cell_organ['cell_id'] = cell_organ['cell_id'].str.upper().str.strip()
cell_organ['organ'] = cell_organ['organ'].str.lower().str.strip()
cell_to_organ = dict(zip(cell_organ['cell_id'], cell_organ['organ']))

scores['cell_organ'] = scores['cell_id'].str.upper().str.strip().map(cell_to_organ)
scores['is_homologous'] = scores['cell_organ'] == scores['target_organ']

# 3. Wilcoxon 秩和检验
print("\n" + "=" * 70)
print("同源 vs 异源 mCS Wilcoxon 秩和检验")
print("=" * 70)

results = []
for organ in sorted(scores['target_organ'].unique()):
    organ_df = scores[scores['target_organ'] == organ].copy()

    homo = organ_df[organ_df['is_homologous'] == True]['mCS'].dropna()
    hetero = organ_df[organ_df['is_homologous'] == False]['mCS'].dropna()

    if len(homo) < 10 or len(hetero) < 10:
        print(f"{organ:12s}: 同源={len(homo)}, 异源={len(hetero)} → 样本不足，跳过")
        continue

    # Wilcoxon 秩和检验
    stat, pval = stats.ranksums(homo, hetero)

    # Cohen's d
    pooled_std = np.sqrt((homo.var() + hetero.var()) / 2)
    cohens_d = (homo.mean() - hetero.mean()) / pooled_std if pooled_std > 0 else 0


    # Cliff's Delta (更稳健的非参数效应量)
    def cliffs_delta(x, y):
        x = np.array(x)
        y = np.array(y)
        nx, ny = len(x), len(y)
        dominance = 0
        for xi in x:
            dominance += np.sum(xi > y) - np.sum(xi < y)
        return dominance / (nx * ny)


    cliff = cliffs_delta(homo, hetero)

    results.append({
        'organ': organ,
        'n_homo': len(homo),
        'n_hetero': len(hetero),
        'mean_homo': homo.mean(),
        'mean_hetero': hetero.mean(),
        'median_homo': homo.median(),
        'median_hetero': hetero.median(),
        'std_homo': homo.std(),
        'std_hetero': hetero.std(),
        'wilcoxon_stat': stat,
        'p_value': pval,
        'cohens_d': cohens_d,
        'cliffs_delta': cliff,
        'significant': pval < 0.05
    })

    sig = "***" if pval < 0.001 else "**" if pval < 0.01 else "*" if pval < 0.05 else "ns"
    print(f"{organ:12s}: homo={len(homo):5d}(μ={homo.mean():.4f}, med={homo.median():.4f}), "
          f"hetero={len(hetero):6d}(μ={hetero.mean():.4f}, med={hetero.median():.4f}), "
          f"d={cohens_d:.3f}, cliff={cliff:.3f}, p={pval:.2e} {sig}")

stat_df = pd.DataFrame(results)
stat_df.to_csv(RESULTS / "homologous_vs_heterologous_stats.csv", index=False)
print(f"\n✓ 保存: {RESULTS / 'homologous_vs_heterologous_stats.csv'}")

# 4. 跨数据集稳健性
print("\n" + "=" * 70)
print("跨数据集稳健性 (Spearman 相关)")
print("=" * 70)

dataset_medians = scores.groupby(['dataset', 'target_organ', 'pert_id', 'pert_iname'])['mCS'].median().reset_index()
pivot = dataset_medians.pivot_table(index=['target_organ', 'pert_id', 'pert_iname'],
                                    columns='dataset', values='mCS').reset_index()

valid = pivot.dropna(subset=['GSE92742', 'GSE70138'])
for organ in sorted(valid['target_organ'].unique()):
    sub = valid[valid['target_organ'] == organ]
    if len(sub) >= 10:
        rho, p = stats.spearmanr(sub['GSE92742'], sub['GSE70138'])
        print(f"{organ:12s}: n={len(sub):5d}, Spearman ρ={rho:.3f}, p={p:.2e}")

print("\n=== 完成 ===")
