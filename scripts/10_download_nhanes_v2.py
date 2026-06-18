import os
import requests
from pathlib import Path
import time

BASE = Path.home() / "project/oadrm/data/raw/nhanes"
BASE.mkdir(parents=True, exist_ok=True)

# 请求头（模拟浏览器）
HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
}

CYCLES = [
    "1999-2000", "2001-2002", "2003-2004", "2005-2006",
    "2007-2008", "2009-2010", "2011-2012", "2013-2014",
    "2015-2016", "2017-2018"
]

def cycle_to_code(cycle):
    parts = cycle.split('-')
    return f"{parts[0][2:]}-{parts[1][2:]}"

def download_file(url, filepath):
    """使用 requests 下载，处理重定向"""
    if filepath.exists():
        # 检查是否是有效 XPT 文件（不是 HTML）
        with open(filepath, 'rb') as f:
            header = f.read(20)
            if b'html' not in header.lower() and b'<!' not in header:
                print(f"  ✓ 已存在且有效: {filepath.name}")
                return True
            else:
                print(f"  ⚠ 存在但无效（HTML），重新下载: {filepath.name}")
    
    try:
        print(f"  下载: {url}")
        response = requests.get(url, headers=HEADERS, allow_redirects=True, timeout=30)
        response.raise_for_status()
        
        # 检查内容类型
        content_type = response.headers.get('content-type', '')
        if 'html' in content_type.lower():
            print(f"  ✗ 返回的是 HTML 页面，不是数据文件")
            return False
        
        # 保存
        with open(filepath, 'wb') as f:
            f.write(response.content)
        
        # 验证
        with open(filepath, 'rb') as f:
            header = f.read(20)
            if b'html' in header.lower() or b'<!' in header:
                print(f"  ✗ 文件内容仍是 HTML")
                return False
        
        print(f"  ✓ 成功: {filepath.name} ({len(response.content)/1024:.1f} KB)")
        return True
        
    except Exception as e:
        print(f"  ✗ 失败: {e}")
        return False

# 需要下载的文件
FILES = {
    "prescription": {
        "filename": "RXQ_RX_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/RXQ_RX_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    "liver": {
        "filename": "BIOPRO_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/BIOPRO_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    "kidney": {
        "filename": "KIDNEY_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/KIDNEY_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    "lung": {
        "filename": "SPX_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/SPX_{cycle_code}.XPT",
        "cycles": ["1999-2000", "2001-2002", "2003-2004", "2005-2006", 
                   "2007-2008", "2009-2010", "2011-2012"]
    },
    "glucose": {
        "filename": "GLU_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/GLU_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    "lipids": {
        "filename": "TCHOL_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/TCHOL_{cycle_code}.XPT",
        "cycles": CYCLES
    },
    "demographics": {
        "filename": "DEMO_{cycle}.XPT",
        "url": "https://wwwn.cdc.gov/Nchs/Nhanes/{cycle}/DEMO_{cycle_code}.XPT",
        "cycles": CYCLES
    }
}

# 执行下载
print("=== NHANES 数据下载（修正版）===")
print(f"保存目录: {BASE}")
print("注意: 如果仍然失败，建议手动从 https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?Cycle=2017-2018 下载\n")

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
        
        time.sleep(0.5)  # 礼貌延迟

print(f"\n=== 下载完成 ===")
print(f"成功: {success_count}, 失败: {fail_count}")

# 清理失败的 HTML 文件
print("\n=== 清理无效文件 ===")
for f in BASE.glob('*.XPT'):
    with open(f, 'rb') as fh:
        header = fh.read(50)
        if b'html' in header.lower() or b'<!' in header or b'DOCTYPE' in header:
            print(f"  删除无效文件: {f.name}")
            f.unlink()
