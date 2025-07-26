import pandas as pd
import os

import scanpy as sc
import pandas as pd
import os

def summarize_cell_gene_counts(input_dir, output_csv_name="samplewise_gene_cell_count.csv"):
    adata_list = {}
    summary = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".h5ad"):
            file_path = os.path.join(input_dir, filename)
            sample_name = os.path.splitext(filename)[0]
            try:
                adata = sc.read_h5ad(file_path)
                adata_list[sample_name] = adata
                summary.append({
                    "sample": sample_name,
                    "n_cells": adata.n_obs,
                    "n_genes": adata.n_vars
                })
            except Exception as e:
                print(f"❌ {sample_name} 읽기 실패: {e}")

    summary_df = pd.DataFrame(summary)
    output_path = os.path.join(input_dir, output_csv_name)
    summary_df.to_csv(output_path, index=False)

    print(f"\n✅ Sample 요약 완료. 결과 저장: {output_path}\n")
    print(summary_df)

    return summary_df

# 📌 여기서 절대경로 지정
input_dir = "/data1/BLA_projection_2/Results/For adata_list(after preprocessing)"
summary_df = summarize_cell_gene_counts(input_dir)


import scanpy as sc
import pandas as pd
import os

# 1. 분석할 단일 adata 파일의 절대 경로
adata_path = "/data1/BLA_projection_2/Results/1_QC/PL_vHPC_1/BLA_projection_PL_vHPC_1_sample_after_preprocessing.h5ad"

# 2. 파일 로딩
adata = sc.read_h5ad(adata_path)

# 3. 샘플 이름 설정 (파일명에서 추출)
sample_name = os.path.splitext(os.path.basename(adata_path))[0]

# 4. 세포 수 / 유전자 수 요약
summary = pd.DataFrame([{
    "sample": sample_name,
    "n_cells": adata.n_obs,
    "n_genes": adata.n_vars
}])

# 5. 출력 및 저장
print(summary)


import scanpy as sc
import os
import pandas as pd
import scrublet as scr 

file_path = "/data1/BLA_projection_2/Results/1_QC/PL_vHPC_3/BLA_projection_sample_PL_vHPC_3.h5ad"

adata = sc.read_h5ad(file_path)

adata = sc.read_h5ad("/data1/BLA_projection_2/Results/1_QC/PL_vHPC_3/BLA_projection_sample_PL_vHPC_3.h5ad")

adata = sc.read_h5ad(file_path)

def run_preprocess(adata_tmp):
    adata_tmp = adata_tmp.copy()

    sample_id = list(adata_tmp.obs['sampleID'].unique())[0]
    print(f'Apply QC for {sample_id}')

    if not pd.api.types.is_string_dtype(adata_tmp.var_names):
        if 'gene_symbols' in adata_tmp.var.columns:
            adata_tmp.var_names = adata_tmp.var['gene_symbols']
        elif 'gene_ids' in adata_tmp.var.columns:
            adata_tmp.var_names = adata_tmp.var['gene_ids']
        adata_tmp.var_names_make_unique()
        adata_tmp.var.index = adata_tmp.var_names

    sc.pp.filter_genes(adata_tmp, min_cells=3)

    adata_tmp = adata_tmp[
        (adata_tmp.obs.n_genes_by_counts >= 500) &
        (adata_tmp.obs.n_genes_by_counts <= 7000) &
        (adata_tmp.obs.total_counts >= 1000) &
        (adata_tmp.obs.pct_counts_mt < 20),
        :
    ]

    try:
        print("📐 Scrublet input shape:", adata_tmp.X.shape)

        scrub = scr.Scrublet(adata_tmp.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False, n_prin_comps=20)

    
        if predicted_doublets is None:
            print("❌ predicted_doublets is None → Scrublet 내부 실패")
        return None
    
        adata_tmp.obs['doublet_scores'] = doublet_scores
        adata_tmp.obs['predicted_doublets'] = predicted_doublets.astype(str)
        adata_tmp = adata_tmp[adata_tmp.obs.doublet_scores < 0.2, :]
    except Exception as e:
        print(f'#Pass {sample_id}: Scrublet failed ({e})')
        return None

    if adata_tmp.n_obs < 500:
        print(f"#Skip {sample_id}: Less than 500 cells after QC")
        return None

    print(f'Preprocessing complete for {sample_id}.')
    return adata_tmp


adata_qc = run_preprocess(adata)

# 4. 결과 저장 (원하면)
if adata_qc is not None:
    save_dir = "/data1/BLA_projection_2/Results/1_QC/PL_vHPC_3"
    filename = "BLA_projection_PL_vHPC_3_sample_after_preprocessing.h5ad"
    save_path = os.path.join(save_dir, filename)
    os.makedirs(save_dir, exist_ok=True)
    adata_qc.write(save_path)
    print(f"✅ Saved to {save_path}")
else:
    print("❌ QC 실패: 저장하지 않음 (None 반환됨)")
    
import scanpy as sc

adata = sc.read_h5ad("/data1/BLA_projection_2/Data/Resources/.v2/after_preprocessing/adata_all_concat.h5ad")

