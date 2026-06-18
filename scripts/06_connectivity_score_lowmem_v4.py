import pandas as pd
import numpy as np
from pathlib import Path
import glob
import h5py
import gc

BASE = Path.home() / "project/oadrm"
SIG_DIR = BASE / "data/processed/organ_signatures"
META_DIR = BASE / "data/processed/lincs_meta"
RAW_LINC = BASE / "data/raw/lincs"
PROC = BASE / "data/processed/drug_scores"
PROC.mkdir(parents=True, exist_ok=True)

BATCH_SIZE = 1000

# ==================== 1. 读取器官签名 ====================
print("=== 读取器官签名 ===")
organs = ['blood', 'brain', 'intestine', 'kidney', 'liver', 'lung', 'skin', 'stomach']
signatures = {}

for organ in organs:
    pattern = f"signature_{organ}_top[0-9]*.csv"
    files = [f for f in glob.glob(str(SIG_DIR / pattern)) if '_protein' not in f and 'summary' not in f]
    if not files:
        continue
    sig = pd.read_csv(files[0])
    sig_dict = {}
    for _, row in sig.iterrows():
        gene = str(row['gene_symbol']).upper()
        if gene and gene != 'NAN':
            if gene not in sig_dict:
                sig_dict[gene] = (float(row['w_norm']), float(row['direction']))
    signatures[organ] = sig_dict
    print(f"  {organ:12s}: {len(sig_dict)} genes")

# ==================== 2. 读取 LINCS 基因映射 ====================
print("\n=== 读取 LINCS 基因信息 ===")
g1 = pd.read_csv(RAW_LINC / "GSE92742_Broad_LINCS_gene_info.txt", sep='\t')
g2 = pd.read_csv(RAW_LINC / "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt", sep='\t')
genes_df = pd.concat([g1, g2]).drop_duplicates(subset='pr_gene_id')
gene_id_to_symbol = dict(zip(
    genes_df['pr_gene_id'].astype(str),
    genes_df['pr_gene_symbol'].str.upper()
))
print(f"LINCS基因数: {len(gene_id_to_symbol)}")

# ==================== 3. 读取元数据 ====================
print("\n=== 读取元数据 ===")
sig_meta = pd.read_csv(META_DIR / "sig_info_trt_cp.csv", low_memory=False)
cell_organ = pd.read_csv(META_DIR / "cellline_to_organ.csv")
valid_cells = set(cell_organ['cell_id'].dropna())
sig_meta = sig_meta[sig_meta['cell_id'].isin(valid_cells)].copy()
sig_meta = sig_meta[['sig_id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose']].copy()

sig_meta['dataset'] = sig_meta['sig_id'].apply(
    lambda x: 'GSE92742' if str(x).startswith('LJP') or str(x).startswith('CPC') else 'GSE70138'
)
print(f"有效记录: {len(sig_meta)}")
print(f"数据集分布:\n{sig_meta['dataset'].value_counts()}")


# ==================== 4. mCS 计算函数 ====================
def compute_mcs_batch(z_matrix, gene_symbols, sig_dict):
    sig_idx = [i for i, g in enumerate(gene_symbols) if g in sig_dict]
    if len(sig_idx) == 0:
        return np.full(z_matrix.shape[0], np.nan)
    sub_z = z_matrix[:, sig_idx]
    weights = np.array([sig_dict[gene_symbols[i]][0] for i in sig_idx])
    directions = np.array([sig_dict[gene_symbols[i]][1] for i in sig_idx])
    scores = np.sum(sub_z * weights * (-directions), axis=1) / len(sig_idx)
    return scores


# ==================== 5. 流式处理（修复版） ====================
def process_gctx_streaming(gctx_path, sig_meta_subset, dataset_name, organ_list):
    print(f"\n{'=' * 50}")
    print(f"处理 {dataset_name}: {gctx_path.name}")
    print(f"{'=' * 50}")

    with h5py.File(gctx_path, 'r') as f:
        # 读取基因符号
        gene_ids = f['0/META/ROW/id'][:]
        gene_symbols = []
        for gid in gene_ids:
            gid_str = gid.decode() if isinstance(gid, bytes) else str(gid)
            gene_symbols.append(gene_id_to_symbol.get(gid_str, 'UNKNOWN'))
        n_genes = len(gene_symbols)

        # 读取样本ID
        sig_ids = f['0/META/COL/id'][:]
        sig_ids = [sid.decode() if isinstance(sid, bytes) else str(sid) for sid in sig_ids]
        n_sigs = len(sig_ids)

        # 检查矩阵形状
        matrix = f['0/DATA/0/matrix']
        print(f"  基因数: {n_genes}, 样本数: {n_sigs}")
        print(f"  矩阵形状: {matrix.shape}")

        # 确定维度布局
        if matrix.shape[0] == n_genes and matrix.shape[1] == n_sigs:
            layout = "genes_x_samples"  # (12328, 473647)
            print(f"  布局: genes x samples")
        elif matrix.shape[0] == n_sigs and matrix.shape[1] == n_genes:
            layout = "samples_x_genes"  # (473647, 12328)
            print(f"  布局: samples x genes")
        else:
            raise ValueError(f"矩阵形状 {matrix.shape} 与预期不符")

        sig_id_to_idx = {sid: i for i, sid in enumerate(sig_ids)}

        subset_meta = sig_meta_subset[sig_meta_subset['sig_id'].isin(sig_id_to_idx)].copy()
        subset_meta = subset_meta.reset_index(drop=True)
        print(f"  匹配样本: {len(subset_meta)}")

        if len(subset_meta) == 0:
            print("  无匹配样本，跳过")
            return

        for organ in organ_list:
            print(f"\n--- 计算 {organ} ---")
            sig_dict = signatures[organ]

            output_file = PROC / f"{dataset_name}_{organ}_mcs.csv"
            n_batches = (len(subset_meta) + BATCH_SIZE - 1) // BATCH_SIZE
            first_batch = True

            for batch_idx in range(n_batches):
                start = batch_idx * BATCH_SIZE
                end = min((batch_idx + 1) * BATCH_SIZE, len(subset_meta))
                batch_meta = subset_meta.iloc[start:end].copy()

                # 获取索引
                indices = [sig_id_to_idx[sid] for sid in batch_meta['sig_id']]

                # 排序并记录恢复顺序
                sorted_indices = sorted(set(indices))  # 去重+排序
                order_map = {idx: i for i, idx in enumerate(sorted_indices)}
                restore_order = [order_map[idx] for idx in indices]

                # 分块读取（避免h5py fancy indexing bug）
                chunk_size = 100
                n_chunks = (len(sorted_indices) + chunk_size - 1) // chunk_size

                if layout == "genes_x_samples":
                    # 矩阵是 (genes, samples)，读取列
                    data = np.zeros((n_genes, len(sorted_indices)), dtype=np.float32)
                    for ci in range(n_chunks):
                        c_start = ci * chunk_size
                        c_end = min((ci + 1) * chunk_size, len(sorted_indices))
                        c_indices = sorted_indices[c_start:c_end]
                        data[:, c_start:c_end] = matrix[:, c_indices]
                    # 恢复顺序并转置为 (samples, genes)
                    data = data[:, restore_order].T
                else:
                    # 矩阵是 (samples, genes)，读取行
                    data = np.zeros((len(sorted_indices), n_genes), dtype=np.float32)
                    for ci in range(n_chunks):
                        c_start = ci * chunk_size
                        c_end = min((ci + 1) * chunk_size, len(sorted_indices))
                        c_indices = sorted_indices[c_start:c_end]
                        data[c_start:c_end, :] = matrix[c_indices, :]
                    # 恢复顺序
                    data = data[restore_order, :]

                # 计算mCS
                scores = compute_mcs_batch(data, gene_symbols, sig_dict)
                batch_meta['mCS'] = scores

                # 保存
                out_cols = ['sig_id', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose', 'mCS']
                out_df = batch_meta[[c for c in out_cols if c in batch_meta.columns]].copy()

                mode = 'w' if first_batch else 'a'
                header = first_batch
                out_df.to_csv(output_file, mode=mode, header=header, index=False)
                first_batch = False

                del data, scores, batch_meta, out_df
                gc.collect()

                if (batch_idx + 1) % 50 == 0 or batch_idx == 0:
                    print(f"  进度: {batch_idx + 1}/{n_batches} ({end}/{len(subset_meta)})")

            print(f"  ✓ 完成: {output_file}")
            gc.collect()


# ==================== 6. 执行 ====================
print("\n" + "=" * 50)
print("开始流式连通性得分计算")
print(f"内存目标: <20GB | Batch size: {BATCH_SIZE}")
print("=" * 50)

gctx1 = RAW_LINC / "GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
meta1 = sig_meta[sig_meta['dataset'] == 'GSE92742']
process_gctx_streaming(gctx1, meta1, "GSE92742", list(signatures.keys()))

gctx2 = RAW_LINC / "GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx"
meta2 = sig_meta[sig_meta['dataset'] == 'GSE70138']
process_gctx_streaming(gctx2, meta2, "GSE70138", list(signatures.keys()))

print("\n=== 全部完成 ===")
print(f"结果目录: {PROC}")
