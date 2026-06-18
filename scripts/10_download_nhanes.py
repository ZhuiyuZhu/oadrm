import os
import urllib.request
from pathlib import Path
import pandas as pd

BASE = Path.home() / "project/oadrm/data/raw/nhanes"
BASE.mkdir(parents=True, exist_ok=True)

# NHANES 周期列表（1999-2018）
CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006",
    "2007-2008", "2009-2010", "2011-2012", "2013-2014",
    "2015-2016", "2017-2018"
]

# 需要下载的文件
FILES = {
    # 用药数据（处方）
    "prescription": {
        "filename": "RXQ_RX_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/RXQ_RX_{cycle_code}.XPT",
        "cycles": CYCLES  # 所有周期都有
    },
    # 肝功能
    "liver": {
        "filename": "BIOPRO_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/BIOPRO_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    # 肾功能
    "kidney": {
        "filename": "KIDNEY_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/KIDNEY_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    # 肺功能（仅到2012）
    "lung": {
        "filename": "SPX_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/SPX_{cycle_code}.XPT",
        "cycles": ["1999-2000", "2001-2002", "2003-2004", "2005-2006",
                   "2007-2008", "2009-2010", "2011-2012"]
    },
    # 血糖
    "glucose": {
        "filename": "GLU_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/GLU_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    # 血脂
    "lipids": {
        "filename": "TCHOL_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/TCHOL_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    # 人口学
    "demographics": {
        "filename": "DEMO_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/DEMO_{cycle_code}.XPT",
        "cycles": CYCLES
    }
}


def cycle_to_code(cycle):
    """1999-2000 -> 99-00"""
    parts = cycle.split('-')
    return f"{parts[0][2:]}-{parts[1][2:]}"


def download_file(url, filepath):
    """下载文件，支持断点续传"""
    if filepath.exists():
        print(f"  已存在: {filepath.name}")
        return True

    try:
        print(f"  下载: {url}")
        urllib.request.urlretrieve(url, str(filepath))
        return True
    except Exception as e:
        print(f"  ✗ 失败: {e}")
        return False


# 执行下载
print("=== NHANES 数据下载 ===")
print(f"保存目录: {BASE}")
print(f"预计下载: {sum(len(v['cycles']) for v in FILES.values())} 个文件")
print("注意: 每个文件约 1-5MB，总下载时间约 10-30 分钟\n")

success_count = 0
fail_count = 0

for category, info in FILES.items():
    print(f"\n--- {category.upper()} ---")
    for cycle in info['cycles']:
        cycle_code = cycle_to_code(cycle)
        filename = info['filename'].format(cycle=cycle)
        url = info['url'].format(cycle=cycle, cycle_code=cycle_code)
        filepath = BASE / filename

        if download_file(url, filepath):
            success_count += 1
        else:
            fail_count += 1

print(f"\n=== 下载完成 ===")
print(f"成功: {success_count}, 失败: {fail_count}")
print(f"文件保存在: {BASE}")
