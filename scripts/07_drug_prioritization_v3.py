import pandas as pd
import numpy as np
from pathlib import Path
import glob

BASE = Path.home() / "project/oadrm"
SCORE_DIR = BASE / "data/processed/drug_scores"
META_DIR = BASE / "data/processed/lincs_meta"
RESULTS = BASE / "results/tables"
RESULTS.mkdir(parents=True, exist_ok=True)

# ==================== 1. 读取 mCS ====================
print("=== 读取 mCS 结果 ===")
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

# ==================== 2. 同源映射 ====================
cell_organ = pd.read_csv(META_DIR / "cellline_to_organ.csv")
cell_organ['cell_id'] = cell_organ['cell_id'].str.upper().str.strip()
cell_organ['organ'] = cell_organ['organ'].str.lower().str.strip()
cell_to_organ = dict(zip(cell_organ['cell_id'], cell_organ['organ']))

scores['cell_organ'] = scores['cell_id'].str.upper().str.strip().map(cell_to_organ)
scores['is_homologous'] = scores['cell_organ'] == scores['target_organ']

# ==================== 3. 核心修正：严格同源过滤 ====================
print("\n=== 核心修正：严格同源过滤 ===")

# 只保留在同源细胞系中有实际数据的药物-器官组合
drug_cell_scores = scores.groupby([
    'pert_id', 'pert_iname', 'cell_id', 'cell_organ',
    'target_organ', 'is_homologous', 'dataset'
])['mCS'].median().reset_index()

# 按药物-器官-数据集聚合同源/异源
homo_scores = drug_cell_scores[drug_cell_scores['is_homologous'] == True].groupby(
    ['pert_id', 'pert_iname', 'target_organ', 'dataset']
)['mCS'].median().reset_index().rename(columns={'mCS': 'mCS_homo'})

hetero_scores = drug_cell_scores[drug_cell_scores['is_homologous'] == False].groupby(
    ['pert_id', 'pert_iname', 'target_organ', 'dataset']
)['mCS'].mean().reset_index().rename(columns={'mCS': 'mCS_hetero'})

# 合并：必须同时有同源和异源数据
organ_scores = homo_scores.merge(
    hetero_scores, on=['pert_id', 'pert_iname', 'target_organ', 'dataset'], how='inner'
)

# 计算 OrganScore（必须同源数据存在）
organ_scores['OrganScore'] = organ_scores['mCS_homo'] - organ_scores['mCS_hetero']

# 跨数据集合并（取中位数）
organ_scores = organ_scores.groupby(['pert_id', 'pert_iname', 'target_organ']).agg({
    'mCS_homo': 'median',
    'mCS_hetero': 'median',
    'OrganScore': 'median'
}).reset_index()

print(f"严格过滤后: {len(organ_scores)} 药物-器官组合")
print(f"各器官组合数:\n{organ_scores['target_organ'].value_counts()}")

# ==================== 4. 四维评分 ====================
print("\n=== 四维评分 ===")

# 4.1 逆转效能
organ_scores['efficacy_z'] = organ_scores.groupby('target_organ')['OrganScore'].transform(
    lambda x: (x - x.mean()) / x.std() if x.std() > 0 else 0
)

# 4.2 特异性
organ_mean = organ_scores.groupby('target_organ')['OrganScore'].transform('mean')
organ_scores['specificity'] = organ_scores['OrganScore'] - organ_mean

# 4.3 稳健性（同源得分稳定性）
organ_scores['robustness'] = np.abs(organ_scores['mCS_homo']) / (
        np.abs(organ_scores['mCS_homo']) + np.abs(organ_scores['mCS_hetero']) + 1e-10
)

# 4.4 综合得分
organ_scores['composite_score'] = (
        0.4 * organ_scores['efficacy_z'] +
        0.3 * organ_scores['specificity'] +
        0.3 * organ_scores['robustness']
)

print(
    f"composite_score: mean={organ_scores['composite_score'].mean():.3f}, max={organ_scores['composite_score'].max():.3f}")

# ==================== 5. 生成候选清单 ====================
print("\n=== 生成候选药物清单 ===")

top_drugs = []
for organ in organ_scores['target_organ'].unique():
    organ_df = organ_scores[organ_scores['target_organ'] == organ].copy()
    organ_df = organ_df.dropna(subset=['composite_score'])
    if len(organ_df) == 0:
        continue
    top20 = organ_df.nlargest(min(20, len(organ_df)), 'composite_score')
    top20['rank'] = range(1, len(top20) + 1)
    top_drugs.append(top20)

if top_drugs:
    top_drugs_df = pd.concat(top_drugs, ignore_index=True)
    print(f"Top 药物: {len(top_drugs_df)} entries ({top_drugs_df['pert_iname'].nunique()} unique)")
    top_drugs_df.to_csv(RESULTS / "organ_drug_top20_v3_filtered.csv", index=False)

    # 广谱 vs 器官特异性
    compound_counts = top_drugs_df.groupby('pert_iname')['target_organ'].nunique().reset_index()
    compound_counts.columns = ['pert_iname', 'n_organs']

    broad = compound_counts[compound_counts['n_organs'] >= 4]['pert_iname'].tolist()
    specific = compound_counts[compound_counts['n_organs'] == 1]['pert_iname'].tolist()

    print(f"\n广谱抗衰老候选 ({len(broad)} 种):")
    for drug in broad[:15]:
        print(f"  {drug}")

    print(f"\n器官特异性候选 ({len(specific)} 种):")
    for drug in specific[:15]:
        print(f"  {drug}")

    # 按器官展示 Top 5
    print("\n=== 每器官 Top 5（严格同源过滤后）===")
    for organ in sorted(top_drugs_df['target_organ'].unique()):
        sub = top_drugs_df[top_drugs_df['target_organ'] == organ].nsmallest(5, 'rank')
        print(f"\n{organ.upper()}:")
        for _, row in sub.iterrows():
            print(f'  #{int(row["rank"])} {row["pert_iname"]:20s} '
                  f'score={row["composite_score"]:.3f} '
                  f'homo={row["mCS_homo"]:.3f} hetero={row["mCS_hetero"]:.3f}')
else:
    print("⚠ 无候选药物")

# 保存完整矩阵
organ_scores.to_csv(RESULTS / "organ_drug_full_matrix_v3.csv", index=False)
print(f"\n✓ 完整矩阵: {RESULTS / 'organ_drug_full_matrix_v3.csv'}")

print("\n=== 完成 ===")
