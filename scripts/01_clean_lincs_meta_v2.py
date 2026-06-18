import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

BASE = Path.home() / "project/oadrm"
RAW = BASE / "data/raw/lincs"
META = BASE / "data/processed/lincs_meta"
META.mkdir(parents=True, exist_ok=True)

# 1. 读取 sig_info
print("=== 读取 GSE92742 sig_info ===")
sig1 = pd.read_csv(RAW / "GSE92742_Broad_LINCS_sig_info.txt", sep="\t", low_memory=False)
print(f"GSE92742: {sig1.shape}")

print("\n=== 读取 GSE70138 sig_info ===")
sig2 = pd.read_csv(RAW / "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt", sep="\t", low_memory=False)
print(f"GSE70138: {sig2.shape}")

# 2. 合并并过滤化合物处理
sig_all = pd.concat([sig1, sig2], ignore_index=True)
print(f"\n合并后: {sig_all.shape}")

sig_cp = sig_all[sig_all['pert_type'] == 'trt_cp'].copy()
print(f"化合物处理 (trt_cp): {sig_cp.shape}")

# 3. 提取关键字段
sig_cp = sig_cp[['sig_id', 'pert_id', 'pert_iname', 'pert_dose', 'pert_dose_unit',
                 'pert_time', 'pert_time_unit', 'cell_id', 'distil_id']].copy()

# 4. 获取所有细胞系
cell_lines = sig_cp['cell_id'].dropna().unique()
print(f"\n涉及细胞系: {len(cell_lines)} 种")
print(f"全部细胞系: {sorted(cell_lines)}")

# 5. 修正版器官映射（基于 LINCS 实际细胞系）
# 注意：LINCS 并非覆盖所有13器官，部分器官无直接对应细胞系
organ_mapping = {
    # 肝脏
    'HEPG2': 'Liver', 'HUH7': 'Liver', 'HEP3B': 'Liver', 'SKHEP1': 'Liver',
    # 肾脏（HEK293T 是胚胎肾细胞系）
    'HEK293': 'Kidney', 'HEK293T': 'Kidney', 'HEKTE': 'Kidney', 'RPTEC': 'Kidney',
    # 心脏（LINCS 中极少心肌细胞系，以下如有则映射）
    'AC16': 'Heart', 'HCM': 'Heart', 'HL1': 'Heart', 'IPSCCM': 'Heart',
    'CM': 'Heart', 'IPSC_CM': 'Heart', 'CARDIOMYOCYTE': 'Heart',
    # 肺
    'A549': 'Lung', 'HCC827': 'Lung', 'H1299': 'Lung', 'H460': 'Lung',
    'NCIH23': 'Lung', 'NCIH358': 'Lung', 'NCIH520': 'Lung', 'NCIH596': 'Lung',
    'NCIH1299': 'Lung', 'NCIH1975': 'Lung', 'NCIH2228': 'Lung', 'NCIH3122': 'Lung',
    'NCIH3255': 'Lung', 'NCIH441': 'Lung', 'NCIH460': 'Lung', 'NCIH647': 'Lung',
    'NCIH1650': 'Lung', 'NCIH1666': 'Lung', 'NCIH1693': 'Lung', 'NCIH1755': 'Lung',
    'NCIH1792': 'Lung', 'NCIH1838': 'Lung', 'NCIH1944': 'Lung', 'NCIH2009': 'Lung',
    'NCIH2030': 'Lung', 'NCIH2073': 'Lung', 'NCIH2122': 'Lung', 'NCIH2170': 'Lung',
    'NCIH2342': 'Lung', 'NCIH2444': 'Lung', 'NCIH292': 'Lung',
    'HCC4006': 'Lung', 'HCC515': 'Lung', 'HCC78': 'Lung', 'HCC15': 'Lung',
    'HCC44': 'Lung', 'HCC1195': 'Lung', 'HCC1359': 'Lung', 'HCC1438': 'Lung',
    'HCC1568': 'Lung', 'HCC1599': 'Lung', 'HCC1833': 'Lung', 'HCC193': 'Lung',
    'HCC2108': 'Lung', 'HCC2279': 'Lung', 'HCC2429': 'Lung', 'HCC2814': 'Lung',
    'HCC2935': 'Lung', 'HCC366': 'Lung', 'HCC4017': 'Lung', 'HCC4032': 'Lung',
    'HCC95': 'Lung', 'HCC1171': 'Lung',
    'CORL23': 'Lung', 'CORL24': 'Lung', 'CORL105': 'Lung', 'CORL279': 'Lung',
    'CALU1': 'Lung', 'CALU3': 'Lung', 'CALU6': 'Lung', 'DV90': 'Lung',
    'EKVX': 'Lung', 'HOP62': 'Lung', 'HOP92': 'Lung', 'NCIH1155': 'Lung',
    'NCIH1299': 'Lung', 'NCIH146': 'Lung', 'NCIH1568': 'Lung', 'NCIH157': 'Lung',
    'NCIH1618': 'Lung', 'NCIH1648': 'Lung', 'NCIH1688': 'Lung', 'NCIH1694': 'Lung',
    'NCIH1703': 'Lung', 'NCIH1755': 'Lung', 'NCIH1770': 'Lung', 'NCIH1792': 'Lung',
    'NCIH1836': 'Lung', 'NCIH1838': 'Lung', 'NCIH1941': 'Lung', 'NCIH1993': 'Lung',
    'NCIH2009': 'Lung', 'NCIH2030': 'Lung', 'NCIH2052': 'Lung', 'NCIH2073': 'Lung',
    'NCIH2085': 'Lung', 'NCIH2087': 'Lung', 'NCIH209': 'Lung', 'NCIH2106': 'Lung',
    'NCIH211': 'Lung', 'NCIH2122': 'Lung', 'NCIH2126': 'Lung', 'NCIH2170': 'Lung',
    'NCIH2172': 'Lung', 'NCIH2195': 'Lung', 'NCIH2196': 'Lung', 'NCIH2227': 'Lung',
    'NCIH2228': 'Lung', 'NCIH226': 'Lung', 'NCIH2286': 'Lung', 'NCIH23': 'Lung',
    'NCIH2342': 'Lung', 'NCIH2405': 'Lung', 'NCIH2444': 'Lung', 'NCIH2468': 'Lung',
    'NCIH250': 'Lung', 'NCIH2522': 'Lung', 'NCIH2541': 'Lung', 'NCIH2596': 'Lung',
    'NCIH2607': 'Lung', 'NCIH2668': 'Lung', 'NCIH2695': 'Lung', 'NCIH2731': 'Lung',
    'NCIH2803': 'Lung', 'NCIH2804': 'Lung', 'NCIH2882': 'Lung', 'NCIH2887': 'Lung',
    'NCIH292': 'Lung', 'NCIH295': 'Lung', 'NCIH2981': 'Lung', 'NCIH3002': 'Lung',
    'NCIH3108': 'Lung', 'NCIH3122': 'Lung', 'NCIH322': 'Lung', 'NCIH3255': 'Lung',
    'NCIH358': 'Lung', 'NCIH378': 'Lung', 'NCIH383': 'Lung', 'NCIH4006': 'Lung',
    'NCIH413': 'Lung', 'NCIH441': 'Lung', 'NCIH446': 'Lung', 'NCIH460': 'Lung',
    'NCIH520': 'Lung', 'NCIH522': 'Lung', 'NCIH526': 'Lung', 'NCIH596': 'Lung',
    'NCIH647': 'Lung', 'NCIH650': 'Lung', 'NCIH661': 'Lung', 'NCIH727': 'Lung',
    'NCIH740': 'Lung', 'NCIH838': 'Lung', 'NCIH889': 'Lung', 'NCIH920': 'Lung',
    'NCIH969': 'Lung', 'NCIH1048': 'Lung', 'NCIH1092': 'Lung', 'NCIH1105': 'Lung',
    'NCIH1155': 'Lung', 'NCIH1184': 'Lung', 'NCIH1238': 'Lung', 'NCIH1293': 'Lung',
    'NCIH1304': 'Lung', 'NCIH1339': 'Lung', 'NCIH1341': 'Lung', 'NCIH1373': 'Lung',
    'NCIH1395': 'Lung', 'NCIH1417': 'Lung', 'NCIH1436': 'Lung', 'NCIH1437': 'Lung',
    'NCIH1450': 'Lung', 'NCIH146': 'Lung', 'NCIH1522': 'Lung', 'NCIH1568': 'Lung',
    'NCIH157': 'Lung', 'NCIH1573': 'Lung', 'NCIH1581': 'Lung', 'NCIH1596': 'Lung',
    'NCIH1618': 'Lung', 'NCIH1648': 'Lung', 'NCIH1650': 'Lung', 'NCIH1666': 'Lung',
    'NCIH1693': 'Lung', 'NCIH1694': 'Lung', 'NCIH1703': 'Lung', 'NCIH1734': 'Lung',
    'NCIH1770': 'Lung', 'NCIH1792': 'Lung', 'NCIH1836': 'Lung', 'NCIH1838': 'Lung',
    'NCIH1869': 'Lung', 'NCIH187': 'Lung', 'NCIH1915': 'Lung', 'NCIH1930': 'Lung',
    'NCIH1941': 'Lung', 'NCIH1944': 'Lung', 'NCIH1963': 'Lung', 'NCIH1975': 'Lung',
    'NCIH1993': 'Lung', 'NCIH2009': 'Lung', 'NCIH2030': 'Lung', 'NCIH2052': 'Lung',
    'NCIH2073': 'Lung', 'NCIH2085': 'Lung', 'NCIH2087': 'Lung', 'NCIH209': 'Lung',
    'NCIH2106': 'Lung', 'NCIH211': 'Lung', 'NCIH2122': 'Lung', 'NCIH2126': 'Lung',
    'NCIH2170': 'Lung', 'NCIH2172': 'Lung', 'NCIH2195': 'Lung', 'NCIH2196': 'Lung',
    'NCIH2227': 'Lung', 'NCIH2228': 'Lung', 'NCIH226': 'Lung', 'NCIH2286': 'Lung',
    'NCIH2405': 'Lung', 'NCIH2444': 'Lung', 'NCIH2468': 'Lung', 'NCIH250': 'Lung',
    'NCIH2522': 'Lung', 'NCIH2541': 'Lung', 'NCIH2596': 'Lung', 'NCIH2607': 'Lung',
    'NCIH2668': 'Lung', 'NCIH2695': 'Lung', 'NCIH2731': 'Lung', 'NCIH2803': 'Lung',
    'NCIH2804': 'Lung', 'NCIH2882': 'Lung', 'NCIH2887': 'Lung', 'NCIH2981': 'Lung',
    'NCIH3002': 'Lung', 'NCIH3108': 'Lung', 'NCIH322': 'Lung', 'NCIH3255': 'Lung',
    'NCIH358': 'Lung', 'NCIH378': 'Lung', 'NCIH383': 'Lung', 'NCIH4006': 'Lung',
    'NCIH413': 'Lung', 'NCIH441': 'Lung', 'NCIH446': 'Lung', 'NCIH460': 'Lung',
    'NCIH520': 'Lung', 'NCIH522': 'Lung', 'NCIH526': 'Lung', 'NCIH596': 'Lung',
    'NCIH647': 'Lung', 'NCIH650': 'Lung', 'NCIH661': 'Lung', 'NCIH727': 'Lung',
    'NCIH740': 'Lung', 'NCIH838': 'Lung', 'NCIH889': 'Lung', 'NCIH920': 'Lung',
    'NCIH969': 'Lung',
    # 脑/神经
    'NPC': 'Brain', 'NEU': 'Brain', 'SHSY5Y': 'Brain', 'SKNSH': 'Brain',
    'U87MG': 'Brain', 'LN229': 'Brain', 'SF268': 'Brain', 'SF295': 'Brain',
    'SF539': 'Brain', 'SNB19': 'Brain', 'SNB75': 'Brain', 'YKG1': 'Brain',
    'SKNMC': 'Brain', 'SKNBE2': 'Brain', 'BE2C': 'Brain', 'CHP134': 'Brain',
    'CHP212': 'Brain', 'COGN415': 'Brain', 'COGN453': 'Brain', 'COGN471': 'Brain',
    'DAOY': 'Brain', 'GIC11': 'Brain', 'GIC20': 'Brain', 'GIC6': 'Brain',
    'GIC8': 'Brain', 'KELLY': 'Brain', 'LAN1': 'Brain', 'NB1': 'Brain',
    'NB1643': 'Brain', 'NB69': 'Brain', 'NLF': 'Brain', 'NMB': 'Brain',
    'SJRH30': 'Brain', 'SJSA1': 'Brain', 'TC32': 'Brain', 'TC71': 'Brain',
    'TR14': 'Brain', 'U251MG': 'Brain', 'UW228': 'Brain', 'YH13': 'Brain',
    # 肌肉
    'HSMM': 'Muscle', 'HSMM2': 'Muscle', 'SKMC': 'Muscle', 'RHAB': 'Muscle',
    'LHCNM2': 'Muscle', 'C25CL48': 'Muscle',
    # 胰腺
    'PANC1': 'Pancreas', 'CAPAN2': 'Pancreas', 'HPAFII': 'Pancreas',
    'BXPC3': 'Pancreas', 'MIAPACA2': 'Pancreas', 'SU8686': 'Pancreas',
    'ASPC1': 'Pancreas', 'CFPAC1': 'Pancreas', 'Capan1': 'Pancreas',
    'PANC0403': 'Pancreas', 'PANC1005': 'Pancreas', 'PANC193': 'Pancreas',
    'PANC210': 'Pancreas', 'PANC3_27': 'Pancreas', 'PANC5_04': 'Pancreas',
    'PANC8_13': 'Pancreas', 'SW1990': 'Pancreas', 'TCCPAN2': 'Pancreas',
    # 乳腺
    'MCF7': 'Breast', 'MCF10A': 'Breast', 'HS578T': 'Breast', 'BT20': 'Breast',
    'MDAMB231': 'Breast', 'MDAMB468': 'Breast', 'T47D': 'Breast', 'SKBR3': 'Breast',
    'BT474': 'Breast', 'BT483': 'Breast', 'BT549': 'Breast', 'CAMA1': 'Breast',
    'HCC1143': 'Breast', 'HCC1187': 'Breast', 'HCC1395': 'Breast', 'HCC1419': 'Breast',
    'HCC1428': 'Breast', 'HCC1500': 'Breast', 'HCC1569': 'Breast', 'HCC1599': 'Breast',
    'HCC1806': 'Breast', 'HCC1937': 'Breast', 'HCC1954': 'Breast', 'HCC202': 'Breast',
    'HCC2157': 'Breast', 'HCC2218': 'Breast', 'HCC38': 'Breast', 'HCC70': 'Breast',
    'MDAMB134VI': 'Breast', 'MDAMB175VII': 'Breast', 'MDAMB361': 'Breast',
    'MDAMB415': 'Breast', 'MDAMB436': 'Breast', 'MDAMB453': 'Breast',
    'UACC812': 'Breast', 'UACC893': 'Breast', 'ZR751': 'Breast', 'ZR7530': 'Breast',
    # 前列腺
    'PC3': 'Prostate', 'LNCAP': 'Prostate', 'VCAP': 'Prostate', 'DU145': 'Prostate',
    '22RV1': 'Prostate', 'NCIH660': 'Prostate', 'LAPC4': 'Prostate',
    # 结肠/肠
    'HCT116': 'Intestine', 'HT29': 'Intestine', 'COLO205': 'Intestine',
    'SW480': 'Intestine', 'SW620': 'Intestine', 'CACO2': 'Intestine',
    'RKO': 'Intestine', 'COLO320': 'Intestine', 'HCC2998': 'Intestine',
    'KM12': 'Intestine', 'DLD1': 'Intestine', 'LOVO': 'Intestine',
    'CL34': 'Intestine', 'COLO201': 'Intestine', 'COLO206F': 'Intestine',
    'COLO320DM': 'Intestine', 'COLO320HSR': 'Intestine', 'COLO678': 'Intestine',
    'COLO741': 'Intestine', 'COLO818': 'Intestine', 'HCC56': 'Intestine',
    'HCT15': 'Intestine', 'HCT8': 'Intestine', 'HT55': 'Intestine',
    'IS1': 'Intestine', 'IS2': 'Intestine', 'IS3': 'Intestine', 'IS4': 'Intestine',
    'LS1034': 'Intestine', 'LS123': 'Intestine', 'LS174T': 'Intestine',
    'LS180': 'Intestine', 'LS411N': 'Intestine', 'LS513': 'Intestine',
    'NCIH508': 'Intestine', 'NCIH716': 'Intestine', 'SNUC1': 'Intestine',
    'SNUC2A': 'Intestine', 'SNUC2B': 'Intestine', 'SNUC5': 'Intestine',
    'SNU61': 'Intestine', 'SW1116': 'Intestine', 'SW1417': 'Intestine',
    'SW1463': 'Intestine', 'SW48': 'Intestine', 'SW837': 'Intestine',
    'T84': 'Intestine', 'WIDR': 'Intestine',
    # 胃
    'AGS': 'Stomach', 'KATOIII': 'Stomach', 'SNU16': 'Stomach', 'MKN45': 'Stomach',
    'MKN1': 'Stomach', 'MKN7': 'Stomach', 'MKN74': 'Stomach', 'NCIH87': 'Stomach',
    'NUGC3': 'Stomach', 'NUGC4': 'Stomach', 'OCUM1': 'Stomach', 'SNU1': 'Stomach',
    'SNU5': 'Stomach', 'SNU620': 'Stomach', 'SNU638': 'Stomach', 'SNU668': 'Stomach',
    # 皮肤
    'A375': 'Skin', 'SKMEL28': 'Skin', 'SKMEL2': 'Skin', 'SKMEL5': 'Skin',
    'SKMEL30': 'Skin', 'MALME3M': 'Skin', 'LOXIMVI': 'Skin', 'WM793': 'Skin',
    'WM88': 'Skin', 'WM115': 'Skin', 'WM2664': 'Skin', 'WM983B': 'Skin',
    'A101D': 'Skin', 'A2058': 'Skin', 'A375P': 'Skin', 'C32': 'Skin',
    'COLO679': 'Skin', 'COLO800': 'Skin', 'COLO829': 'Skin', 'COLO849': 'Skin',
    'HS294T': 'Skin', 'HS695T': 'Skin', 'HS839T': 'Skin', 'HS936T': 'Skin',
    'IGR39': 'Skin', 'IGR37': 'Skin', 'MEWO': 'Skin', 'RPMI7951': 'Skin',
    'SKMEL1': 'Skin', 'SKMEL24': 'Skin', 'SKMEL3': 'Skin', 'SKMEL31': 'Skin',
    # 卵巢
    'OV7': 'Ovary', 'IGROV1': 'Ovary', 'OVCAR3': 'Ovary', 'OVCAR4': 'Ovary',
    'OVCAR5': 'Ovary', 'OVCAR8': 'Ovary', 'SKOV3': 'Ovary', 'TOV21G': 'Ovary',
    'A2780': 'Ovary', 'CAOV3': 'Ovary', 'COV318': 'Ovary', 'COV362': 'Ovary',
    'COV413': 'Ovary', 'COV434': 'Ovary', 'COV504': 'Ovary', 'COV644': 'Ovary',
    'EFO21': 'Ovary', 'EFO27': 'Ovary', 'FUOV1': 'Ovary', 'HEY': 'Ovary',
    'JHOS2': 'Ovary', 'JHOS4': 'Ovary', 'KURAMOCHI': 'Ovary', 'NIHOVCAR3': 'Ovary',
    'OAW28': 'Ovary', 'OAW42': 'Ovary', 'OV90': 'Ovary', 'OVCA420': 'Ovary',
    'OVCA432': 'Ovary', 'OVCA433': 'Ovary', 'OVCA429': 'Ovary', 'OVCA429': 'Ovary',
    'OVISE': 'Ovary', 'OVMANA': 'Ovary', 'OVTOKO': 'Ovary', 'PEO1': 'Ovary',
    'PEO14': 'Ovary', 'PEO16': 'Ovary', 'PEO23': 'Ovary', 'PEO4': 'Ovary',
    'PEO6': 'Ovary', 'RMUGS': 'Ovary', 'RMUGL': 'Ovary', 'SHIN3': 'Ovary',
    'SNU8': 'Ovary', 'SNU119': 'Ovary', 'SNU251': 'Ovary', 'SNU840': 'Ovary',
    'TYKNU': 'Ovary', 'UWB1B289': 'Ovary',
    # 宫颈
    'HELA': 'Cervix', 'C33A': 'Cervix', 'SIHA': 'Cervix', 'CASKI': 'Cervix',
    'C4II': 'Cervix', 'C4I': 'Cervix', 'DOVK': 'Cervix', 'HCC94': 'Cervix',
    'HT3': 'Cervix', 'ME180': 'Cervix', 'MS751': 'Cervix', 'SW756': 'Cervix',
    # 骨骼
    'U2OS': 'Bone', 'SAOS2': 'Bone', 'MG63': 'Bone', 'HOS': 'Bone',
    '143B': 'Bone', 'CAL72': 'Bone', 'G292': 'Bone', 'HSOS1': 'Bone',
    'HUO3N1': 'Bone', 'HUO9': 'Bone', 'KPD': 'Bone', 'MNNG': 'Bone',
    'NY': 'Bone', 'OS152': 'Bone', 'OS187': 'Bone', 'OS252': 'Bone',
    'OS296': 'Bone', 'OS332': 'Bone', 'OS33': 'Bone', 'OS349': 'Bone',
    'OS36': 'Bone', 'OS384': 'Bone', 'OS391': 'Bone', 'OS402': 'Bone',
    'OS447': 'Bone', 'OS452': 'Bone', 'OS482': 'Bone', 'OS518': 'Bone',
    'OS556': 'Bone', 'OS625': 'Bone', 'OS649': 'Bone', 'OS687': 'Bone',
    'OS74': 'Bone', 'OS832': 'Bone', 'OS9': 'Bone', 'SJSA1': 'Bone',
    'TE85': 'Bone', 'U2OSR1': 'Bone', 'U2OSR2': 'Bone', 'U2OSR3': 'Bone',
    'U2OSR4': 'Bone', 'ZOS': 'Bone',
    # 血液/免疫
    'THP1': 'Blood', 'HL60': 'Blood', 'K562': 'Blood', 'MV411': 'Blood',
    'MOLT4': 'Blood', 'JURKAT': 'Blood', 'CEMC7': 'Blood', 'RPMI8226': 'Blood',
    'U937': 'Blood', 'AMO1': 'Blood', 'MM1S': 'Blood', 'OCIAML2': 'Blood',
    'OCIAML3': 'Blood', 'PL21': 'Blood', 'KG1': 'Blood', 'KASUMI2': 'Blood',
    'KASUMI3': 'Blood', 'MOLM13': 'Blood', 'MOLM14': 'Blood', 'NB4': 'Blood',
    'OCIAML5': 'Blood', 'REH': 'Blood', 'RS411': 'Blood', 'SHI1': 'Blood',
    'SKM1': 'Blood', 'TALL1': 'Blood', 'TF1': 'Blood', 'UMCAL1': 'Blood',
    'BDCM': 'Blood', 'CD34': 'Blood', 'EM2': 'Blood', 'EM3': 'Blood',
    'F36P': 'Blood', 'HEL': 'Blood', 'KU812': 'Blood', 'LAMA84': 'Blood',
    'ME1': 'Blood', 'MUTZ1': 'Blood', 'MUTZ8': 'Blood', 'NOMO1': 'Blood',
    'OCIAML1': 'Blood', 'OCIAML6': 'Blood', 'SKNO1': 'Blood', 'UT7': 'Blood',
    # 甲状腺
    '8505C': 'Thyroid', 'BHT101': 'Thyroid', 'CAL62': 'Thyroid', 'C643': 'Thyroid',
    'KTC1': 'Thyroid', 'SW1736': 'Thyroid', 'TPC1': 'Thyroid',
    'BCPAP': 'Thyroid', 'FTC133': 'Thyroid', 'FTC238': 'Thyroid', 'FTC238': 'Thyroid',
    'KTC2': 'Thyroid', 'RO82W1': 'Thyroid', 'TT': 'Thyroid', 'WRO': 'Thyroid',
    # 肾上腺
    'SW13': 'Adrenal', 'NCIH295R': 'Adrenal',
    # 脾脏（LINCS中无直接脾脏细胞系，可用血液细胞系近似或标记为NA）
    # 膀胱
    'RT4': 'Bladder', 'SW780': 'Bladder', 'T24': 'Bladder', 'UMUC3': 'Bladder',
    '5637': 'Bladder', 'HT1197': 'Bladder', 'HT1376': 'Bladder', 'J82': 'Bladder',
    'SCABER': 'Bladder', 'KU1919': 'Bladder', 'TCCSUP': 'Bladder',
    # 子宫内膜
    'HEC1A': 'Endometrium', 'HEC1B': 'Endometrium', 'HEC50B': 'Endometrium',
    'HEC59': 'Endometrium', 'HEC108': 'Endometrium', 'HEC151': 'Endometrium',
    'HEC265': 'Endometrium', 'HEC6': 'Endometrium', 'ISHIKAWA': 'Endometrium',
    'KLE': 'Endometrium', 'RL952': 'Endometrium', 'SPEC2': 'Endometrium',
    'AN3CA': 'Endometrium', 'ARK1': 'Endometrium', 'EMTK': 'Endometrium',
    'ESS1': 'Endometrium', 'HEC116': 'Endometrium', 'HEC251': 'Endometrium',
    'HEC253': 'Endometrium', 'HEC265B': 'Endometrium', 'HEC88NU': 'Endometrium',
    'HEIKA': 'Endometrium', 'SNGM': 'Endometrium', 'SNU685': 'Endometrium',
    # 食管
    'KYSE30': 'Esophagus', 'KYSE70': 'Esophagus', 'KYSE150': 'Esophagus',
    'KYSE180': 'Esophagus', 'KYSE270': 'Esophagus', 'KYSE410': 'Esophagus',
    'KYSE450': 'Esophagus', 'KYSE510': 'Esophagus', 'TE1': 'Esophagus',
    'TE10': 'Esophagus', 'TE11': 'Esophagus', 'TE12': 'Esophagus',
    'TE15': 'Esophagus', 'TE2': 'Esophagus', 'TE4': 'Esophagus', 'TE5': 'Esophagus',
    'TE6': 'Esophagus', 'TE8': 'Esophagus', 'TE9': 'Esophagus',
    # 头颈
    'CAL27': 'HeadNeck', 'CAL33': 'HeadNeck', 'FADU': 'HeadNeck', 'SCC9': 'HeadNeck',
    'SCC15': 'HeadNeck', 'SCC25': 'HeadNeck', 'SCC4': 'HeadNeck', 'SCC7': 'HeadNeck',
    'TU167': 'HeadNeck', 'UPCI_SCC090': 'HeadNeck', 'UPCI_SCC152': 'HeadNeck',
    'UPCI_SCC154': 'HeadNeck', 'UPCI_SCC172': 'HeadNeck', 'UPCI_SCC200': 'HeadNeck',
    'UPCI_SCC204': 'HeadNeck', 'UPCI_SCC217': 'HeadNeck', 'UPCI_SCC243': 'HeadNeck',
    'UPCI_SCC266': 'HeadNeck', 'UPCI_SCC299': 'HeadNeck', 'UPCI_SCC374': 'HeadNeck',
    'UPCI_SCC399': 'HeadNeck', 'UPCI_SCC452': 'HeadNeck', 'UPCI_SCC497': 'HeadNeck',
    'UPCI_SCC523': 'HeadNeck', 'UPCI_SCC525': 'HeadNeck', 'UPCI_SCC553': 'HeadNeck',
    'UPCI_SCC584': 'HeadNeck', 'UPCI_SCC686': 'HeadNeck', 'UPCI_SCC692': 'HeadNeck',
    'UPCI_SCC740': 'HeadNeck', 'UPCI_SCC745': 'HeadNeck', 'UPCI_SCC782': 'HeadNeck',
    'UPCI_SCC798': 'HeadNeck', 'UPCI_SCC804': 'HeadNeck', 'UPCI_SCC821': 'HeadNeck',
    'UPCI_SCC842': 'HeadNeck', 'UPCI_SCC896': 'HeadNeck', 'UPCI_SCC921': 'HeadNeck',
    'UPCI_SCC939': 'HeadNeck', 'UPCI_SCC942': 'HeadNeck', 'UPCI_SCC953': 'HeadNeck',
    'UPCI_SCC966': 'HeadNeck', 'UPCI_SCC975': 'HeadNeck', 'UPCI_SCC978': 'HeadNeck',
    'UPCI_SCC993': 'HeadNeck', 'UPCI_SCC996': 'HeadNeck', 'UPCI_SCC999': 'HeadNeck',
    'BHY': 'HeadNeck', 'CAL33': 'HeadNeck', 'HSC2': 'HeadNeck', 'HSC3': 'HeadNeck',
    'HSC4': 'HeadNeck', 'PECAPJ34CL': 'HeadNeck', 'PECAPJ41CL': 'HeadNeck',
    'PECAPJ49CL': 'HeadNeck', 'SCC9': 'HeadNeck', 'YD8': 'HeadNeck', 'YD10B': 'HeadNeck',
    'YD38': 'HeadNeck', 'SCC25': 'HeadNeck',
    # 眼
    'Y79': 'Eye', 'WERIRB1': 'Eye',
    # 胎盘
    'JEG3': 'Placenta', 'JAR': 'Placenta', 'BEWO': 'Placenta',
    # 睾丸
    'NT2D1': 'Testis', 'TCAM2': 'Testis', 'NTERA2': 'Testis',
    # 其他常见
    'HA1E': 'Kidney',  # HA1E 是胚胎肾来源
    'HEK293': 'Kidney', 'HEK293T': 'Kidney',
    'CL33': 'Lung', 'CL40': 'Lung',
    'WSUDLCL2': 'Blood', 'WSUDLCL3': 'Blood', 'WSUDLCL4': 'Blood',
    'WSUDLCL5': 'Blood', 'WSUDLCL6': 'Blood', 'WSUDLCL7': 'Blood',
    'WSUDLCL8': 'Blood', 'WSUDLCL9': 'Blood', 'WSUDLCL10': 'Blood',
}

# 6. 应用映射
sig_cp['organ'] = sig_cp['cell_id'].map(organ_mapping)

# 7. 诊断输出
print(f"\n=== 映射后器官分布 ===")
organ_counts = sig_cp['organ'].value_counts(dropna=False)
print(organ_counts)

print(f"\n=== 未映射细胞系 ({sig_cp['organ'].isna().sum()} 条记录) ===")
unmapped = sig_cp[sig_cp['organ'].isna()]['cell_id'].value_counts()
print(unmapped.head(20))

# 8. 保存为 CSV（修复 parquet 问题）
sig_cp.to_csv(META / "sig_info_trt_cp.csv", index=False)
print(f"\n✓ 保存至: {META / 'sig_info_trt_cp.csv'}")

# 9. 生成细胞系-器官映射表
cell_organ = sig_cp[['cell_id', 'organ']].dropna().drop_duplicates()
cell_organ.to_csv(META / "cellline_to_organ.csv", index=False)
print(f"✓ 细胞系-器官映射: {META / 'cellline_to_organ.csv'} ({len(cell_organ)} 条)")

# 10. 生成化合物列表
compounds = sig_cp[['pert_id', 'pert_iname']].drop_duplicates()
compounds.to_csv(META / "compound_list.csv", index=False)
print(f"✓ 化合物列表: {META / 'compound_list.csv'} ({len(compounds)} 种)")

# 11. 生成13器官覆盖报告
target_organs = ['Heart', 'Liver', 'Kidney', 'Brain', 'Lung', 'Muscle', 
                   'Pancreas', 'Stomach', 'Intestine', 'Spleen', 
                   'Adrenal', 'Thyroid', 'Bone']
print(f"\n=== 13器官覆盖报告 ===")
for organ in target_organs:
    count = organ_counts.get(organ, 0)
    status = "✓" if count > 0 else "✗ 缺失"
    print(f"  {organ:12s}: {count:8d} 条记录  {status}")
