import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import stats
import statsmodels.formula.api as smf
import warnings

warnings.filterwarnings('ignore')

BASE = Path.home() / "project/oadrm"
NHANES_DIR = BASE / "data/raw/nhanes/manual_download"
RESULTS = BASE / "results/tables"
FIGURES = BASE / "results/figures"
FIGURES.mkdir(parents=True, exist_ok=True)

print("=" * 60)
print("NHANES 真实世界验证 (2017-2018 Cycle)")
print("=" * 60)

# 1. 读取数据
print("\n=== 1. 读取 NHANES 数据 ===")
demo = pd.read_sas(NHANES_DIR / "DEMO_J.XPT", format='xport')
rx = pd.read_sas(NHANES_DIR / "RXQ_RX_J.XPT", format='xport')
bio = pd.read_sas(NHANES_DIR / "BIOPRO_J.XPT", format='xport')
alb = pd.read_sas(NHANES_DIR / "ALB_CR_J.XPT", format='xport')
print(f"DEMO: {demo.shape}, RXQ_RX: {rx.shape}, BIOPRO: {bio.shape}, ALB_CR: {alb.shape}")

# 2. 药物匹配
print("\n=== 2. 药物匹配 ===")
CANDIDATES = {
    'metformin': {'organ': 'liver', 'type': 'biguanide', 'lincs_predicted': False},
    'atorvastatin': {'organ': 'liver', 'type': 'statin', 'lincs_predicted': False},
    'lisinopril': {'organ': 'kidney', 'type': 'ACE inhibitor', 'lincs_predicted': False},
    'enalapril': {'organ': 'kidney', 'type': 'ACE inhibitor', 'lincs_predicted': False},
    'losartan': {'organ': 'kidney', 'type': 'ARB', 'lincs_predicted': False},
    'simvastatin': {'organ': 'liver', 'type': 'statin', 'lincs_predicted': False},
    'quinapril': {'organ': 'brain', 'type': 'ACE inhibitor', 'lincs_predicted': True},
}

rx['drug_name'] = rx['RXDDRUG'].astype(str).str.lower().str.strip()

matched_drugs = {}
for drug_key, info in CANDIDATES.items():
    mask = rx['drug_name'].str.contains(drug_key, na=False, case=False, regex=False)
    users = rx.loc[mask, 'SEQN'].unique()
    matched_drugs[drug_key] = {
        'users': set(users), 'n': len(users), 'organ': info['organ'],
        'type': info['type'], 'lincs_predicted': info['lincs_predicted']
    }
    print(f"  {drug_key:20s}: {len(users):4d} users ({info['organ']}) {'[LINCS]' if info['lincs_predicted'] else ''}")

# 3. 构建数据集
print("\n=== 3. 构建分析数据集 ===")
df = demo[['SEQN', 'RIAGENDR', 'RIDAGEYR']].copy()
df.columns = ['SEQN', 'sex', 'age']

# 肝功能
bio_cols = ['SEQN', 'LBXSATSI', 'LBXSASSI', 'LBXSAPSI']
bio_cols = [c for c in bio_cols if c in bio.columns]
if bio_cols:
    df = df.merge(bio[bio_cols], on='SEQN', how='left')
    df = df.rename(columns={'LBXSATSI': 'ALT', 'LBXSASSI': 'AST', 'LBXSAPSI': 'ALP'})

# 肾功能
if 'URXUMA' in alb.columns and 'URXUCR' in alb.columns:
    df = df.merge(alb[['SEQN', 'URXUMA', 'URXUCR']], on='SEQN', how='left')
    df['UACR'] = (df['URXUMA'] / df['URXUCR'] * 100).replace([np.inf, -np.inf], np.nan)

# 药物标志
for drug_key, info in matched_drugs.items():
    df[f'use_{drug_key}'] = df['SEQN'].isin(info['users']).astype(int)

print(f"分析数据集: {df.shape}")
print(f"列: {df.columns.tolist()}")

# 4. 关联分析
print("\n=== 4. 药物-器官功能关联分析 ===")

# 只分析有实际功能指标的器官
ORGAN_MARKERS = {
    'liver': ['ALT', 'AST', 'ALP'],
    'kidney': ['UACR'],
}

results = []

for drug_key, info in matched_drugs.items():
    if info['n'] < 30:
        continue

    drug_col = f'use_{drug_key}'
    organ = info['organ']
    markers = [m for m in ORGAN_MARKERS.get(organ, []) if m in df.columns]

    if not markers:
        print(f"  {drug_key}: 跳过 (无直接{organ}功能指标)")
        continue

    print(f"\n--- {drug_key} ({organ}, n={info['n']}) ---")

    for marker in markers:
        # 避免列名冲突：marker 不能是 'age' 或 'sex'
        if marker in ['age', 'sex']:
            continue

        sub = df[[drug_col, marker, 'age', 'sex']].dropna()
        if len(sub) < 30:
            continue

        users = sub[sub[drug_col] == 1][marker]
        nonusers = sub[sub[drug_col] == 0][marker]

        if len(users) < 10 or len(nonusers) < 10:
            continue

        print(f"  {marker}: users={users.mean():.2f}±{users.std():.2f}, "
              f"nonusers={nonusers.mean():.2f}±{nonusers.std():.2f}")

        stat, pval = stats.ranksums(users, nonusers)

        try:
            model = smf.ols(f'{marker} ~ {drug_col} + age + sex', data=sub).fit()
            beta = model.params[drug_col]
            se = model.bse[drug_col]
            p_reg = model.pvalues[drug_col]
            ci_low, ci_high = model.conf_int().loc[drug_col]
        except:
            beta = se = p_reg = ci_low = ci_high = np.nan

        results.append({
            'drug': drug_key, 'organ': organ, 'marker': marker,
            'n_users': len(users), 'n_nonusers': len(nonusers),
            'mean_users': users.mean(), 'mean_nonusers': nonusers.mean(),
            'wilcoxon_p': pval, 'beta': beta, 'se': se,
            'p_regression': p_reg, 'ci_low': ci_low, 'ci_high': ci_high,
            'lincs_predicted': info['lincs_predicted']
        })

if results:
    result_df = pd.DataFrame(results)
    result_df.to_csv(RESULTS / "nhanes_validation_results.csv", index=False)
    print(f"\n✓ 结果保存: {RESULTS / 'nhanes_validation_results.csv'}")
    print(result_df.to_string())
else:
    print("\n⚠ 无有效分析结果")

# 5. 可视化
print("\n=== 5. 可视化 ===")

if results:
    result_df = pd.DataFrame(results)

    fig, ax = plt.subplots(figsize=(10, max(6, len(result_df) * 0.5)))
    y_pos = np.arange(len(result_df))

    colors = []
    for _, row in result_df.iterrows():
        if row['lincs_predicted']:
            colors.append('blue' if row['beta'] < 0 and row['p_regression'] < 0.05 else 'lightblue')
        else:
            colors.append('green' if row['beta'] < 0 and row['p_regression'] < 0.05 else
                          'red' if row['beta'] > 0 and row['p_regression'] < 0.05 else 'gray')

    ax.errorbar(result_df['beta'], y_pos,
                xerr=[result_df['beta'] - result_df['ci_low'],
                      result_df['ci_high'] - result_df['beta']],
                fmt='o', color='black', ecolor=colors, capsize=3, markersize=8)

    ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    ax.set_yticks(y_pos)
    labels = [f"{r['drug']}{' [LINCS]' if r['lincs_predicted'] else ''}\n({r['marker']})"
              for _, r in result_df.iterrows()]
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel('Effect Size (β, drug users vs non-users)', fontsize=11)
    ax.set_title('NHANES Validation: Drug-Organ Function Associations\n'
                 '(Negative β = protective; corrected for age, sex)', fontsize=12)

    for i, (_, row) in enumerate(result_df.iterrows()):
        sig = '***' if row['p_regression'] < 0.001 else '**' if row['p_regression'] < 0.01 else '*' if row[
                                                                                                           'p_regression'] < 0.05 else 'ns'
        ax.text(row['ci_high'] + 0.05, i, sig, va='center', fontsize=9)

    from matplotlib.patches import Patch

    legend_elements = [
        Patch(facecolor='green', label='Known protective (p<0.05)'),
        Patch(facecolor='red', label='Harmful (p<0.05)'),
        Patch(facecolor='blue', label='LINCS predicted (p<0.05)'),
        Patch(facecolor='gray', label='Not significant')
    ]
    ax.legend(handles=legend_elements, loc='lower right', fontsize=8)

    plt.tight_layout()
    plt.savefig(FIGURES / "nhanes_forest_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ 森林图: {FIGURES / 'nhanes_forest_plot.png'}")

print("\n=== 完成 ===")
