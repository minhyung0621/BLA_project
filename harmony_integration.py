import scanpy as sc
import anndata
import os
import glob

# 🔹 루트 디렉토리 설정 (하위 폴더 포함)
data_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_preprocessing"

# 🔹 하위 폴더 포함 .h5ad 파일 전체 검색
adata_paths = sorted(glob.glob(os.path.join(data_dir, "**", "*.h5ad"), recursive=True))

# 🔹 개별 파일 로드
adata_list = [sc.read_h5ad(p) for p in adata_paths]

# 🔹 샘플 이름 추출 (폴더/파일명 기반으로 생성)
sample_keys = [os.path.splitext(os.path.basename(p))[0] for p in adata_paths]

# 🔹 Concat 수행
adata_all = anndata.concat(adata_list, join="outer", label="sample", keys=sample_keys)

# 🔹 저장 경로 설정 및 저장
save_path = os.path.join(data_dir, "adata_all_concat.h5ad")
adata_all.write(save_path)

print(f"✅ Concat 완료 및 저장 경로: {save_path}")

import scanpy as sc

adata_all = sc.read_h5ad("/data1/BLA_projection_2/Data/Resources/.v2/after_preprocessing/adata_all_concat.h5ad")
# Normalize → log1p 적용
sc.pp.normalize_total(adata_all, target_sum=1e4)
sc.pp.log1p(adata_all)

# 미토콘드리아 유전자 비율 계산
adata_all.var['mt'] = adata_all.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_all, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata_all, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], groupby='sample', multi_panel=True)


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

# sample 이름 간단히 정제 (이미 했으면 생략 가능)
adata_all.obs['sample'] = adata_all.obs['sample'].str.replace("BLA_projection_sample_", "", regex=False)

# 필요한 데이터프레임 추출
qc_df = adata_all.obs[['sample', 'n_genes_by_counts', 'pct_counts_mt']].copy()


# 미토콘드리아 유전자 비율 계산
adata_all.var['mt'] = adata_all.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_all, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# QC metric 시각화
sc.pl.violin(adata_all, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], groupby='sample', multi_panel=True)


# PCA 계산 및 layer 저장
sc.pp.highly_variable_genes(adata_all, flavor='seurat', n_top_genes=2000, subset=True)
sc.pp.scale(adata_all, max_value=10)
sc.tl.pca(adata_all, svd_solver='arpack')

# UMAP 계산
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['sample', 'total_counts', 'pct_counts_mt'])


%pip install harmonypy

import harmonypy as hm

# Harmony 적용 전 neighbors 재계산
sc.pp.pca(adata_all, svd_solver='arpack')
ho = hm.run_harmony(adata_all.obsm['X_pca'], adata_all.obs, 'sample')
adata_all.obsm['X_pca_harmony'] = ho.Z_corr.T

# UMAP (Harmony)
sc.pp.neighbors(adata_all, use_rep='X_pca_harmony')
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=['sample', 'total_counts', 'pct_counts_mt'], title='Harmony UMAP')


from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt

# 🔹 다양한 resolution에 대해 실루엣 점수 계산
resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
scores = {}

for res in resolutions:
    # Leiden clustering 수행 (Harmony latent 기반)
    sc.tl.leiden(adata_all, resolution=res, key_added=f'leiden_{res}')
    
    # 실루엣 점수 계산
    labels = adata_all.obs[f'leiden_{res}']
    score = silhouette_score(adata_all.obsm['X_pca_harmony'], labels)
    scores[res] = score
    print(f"🔹 resolution={res:.1f} → silhouette_score={score:.4f}")

# 🔹 시각화
plt.figure(figsize=(6,4))
plt.plot(list(scores.keys()), list(scores.values()), marker='o')
plt.xlabel("Leiden Resolution")
plt.ylabel("Silhouette Score")
plt.title("Resolution vs Silhouette Score")
plt.grid(True)
plt.tight_layout()
plt.show()

# 🔹 최적 resolution으로 Leiden clustering 수행
sc.tl.leiden(adata_all, resolution=0.2, key_added='leiden_0.2')

# 🔹 UMAP 시각화 (클러스터별로 색상 구분)
sc.pl.umap(
    adata_all,
    color='leiden_0.2',
    palette='tab20',  # 색상 팔레트 설정 (최대 20개 클러스터까지)
    legend_loc='on data',  # UMAP 위에 label 표시
    title='Leiden Clustering (resolution=0.2)'
)

# 🔹 클러스터별 marker gene 탐색
sc.tl.rank_genes_groups(
    adata_all,
    groupby='leiden_0.2',
    method='wilcoxon',  # 또는 't-test', 'logreg'
    key_added='rank_genes_leiden_0.2'
)

# 🔹 상위 marker 유전자 시각화
sc.pl.rank_genes_groups(
    adata_all,
    key='rank_genes_leiden_0.2',
    n_genes=10,  # 클러스터당 상위 10개
    sharey=False
)

%pip show celltypist

import celltypist
from celltypist.models import Model

from celltypist.models import Model, download_models

import scanpy as sc

import numpy as np

# 총 카운트를 10,000으로 정규화
sc.pp.normalize_total(adata_all, target_sum=1e4)

# log1p 변환
sc.pp.log1p(adata_all)

adata_all.X = np.nan_to_num(adata_all.X)


# 🔹 1. 모델 다운로드
model_path = download_models(model='Mouse_Whole_Brain.pkl')  # 이름에 공백 ❌, 언더스코어 ✅

# 🔹 2. 모델 로드
model = Model.load('Mouse_Whole_Brain.pkl')

# 🔹 3. CellTypist 예측 실행
results = celltypist.annotate(adata_all, model=model, majority_voting=True)


# 1. 클러스터 단위로 celltypist 실행
results = celltypist.annotate(adata_all, model=model, majority_voting=True)

# 2. 클러스터 단위 결과 → 딕셔너리로 저장 (cluster ID → cell type)
cluster_labels = results.predicted_labels
cluster_to_celltype = dict(zip(cluster_labels.index.astype(str), cluster_labels.values))

# 3. 클러스터 정보를 기반으로 각 cell에 해당 cell type 할당
adata_all.obs['celltypist'] = adata_all.obs['leiden_0.2'].map(cluster_to_celltype)



