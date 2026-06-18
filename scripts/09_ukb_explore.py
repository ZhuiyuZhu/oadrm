import pandas as pd
import numpy as np
from pathlib import Path
import glob

UKB_DIR = Path("/storageB/wangjiahao/data/UKB_phenotype")
BASE = Path.home() / "project/oadrm"
RESULTS = BASE / "results/tables"

print("=== UKB 数据目录结构 ===")
for item in sorted(UKB_DIR.iterdir())[:20]:
    if item.is_dir():
        print(f"  [DIR]  {item.name}")
    elif item.is_file():
        size = item.stat().st_size / 1024 / 1024
        print(f"  [FILE] {item.name} ({size:.1f} MB)")

# 查找用药相关字段
print("\n=== 搜索用药/处方相关文件 ===")
med_files = list(UKB_DIR.glob("*med*")) + list(UKB_DIR.glob("*drug*")) + \
            list(UKB_DIR.glob("*prescript*")) + list(UKB_DIR.glob("*20003*")) + \
            list(UKB_DIR.glob("*6153*")) + list(UKB_DIR.glob("*ATC*"))
for f in med_files[:10]:
    print(f"  {f.name}")

# 查找器官功能指标
print("\n=== 搜索器官功能指标文件 ===")
organ_files = list(UKB_DIR.glob("*eGFR*")) + list(UKB_DIR.glob("*creatinine*")) + \
              list(UKB_DIR.glob("*ALT*")) + list(UKB_DIR.glob("*AST*")) + \
              list(UKB_DIR.glob("*FEV*")) + list(UKB_DIR.glob("*NTproBNP*")) + \
              list(UKB_DIR.glob("*blood*")) + list(UKB_DIR.glob("*liver*")) + \
              list(UKB_DIR.glob("*kidney*")) + list(UKB_DIR.glob("*lung*"))
for f in organ_files[:10]:
    print(f"  {f.name}")

# 如果有CSV/TSV文件，查看列名
print("\n=== 查看样本文件结构 ===")
sample_files = list(UKB_DIR.glob("*.csv")) + list(UKB_DIR.glob("*.tsv")) + list(UKB_DIR.glob("*.txt"))
if sample_files:
    f = sample_files[0]
    print(f"读取: {f.name}")
    try:
        df = pd.read_csv(f, nrows=5)
        print(f"Shape: {df.shape}")
        print(f"Columns: {df.columns.tolist()[:20]}")
    except Exception as e:
        print(f"读取失败: {e}")
