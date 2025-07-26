import scanpy as sc
import scvi

adata = sc.read_h5ad("/data1/BLA_projection_2/Data/Resources/.v2/after_preprocessing/adata_all_concat.h5ad")

adata.layers["counts"] = adata.X.copy()


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.layers["logcounts"] = adata.X.copy()

sc.pp.highly_variable_genes(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
adata

sc.pl.umap(adata, color=["sampleID"], wspace=1)
batch_key = "sampleID"
sc.pl.umap(adata, color=[batch_key], wspace=1)

# 미토콘드리아 유전자 비율 계산
adata.var['mt'] = adata.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# QC metric 시각화
sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt'], groupby='sample', multi_panel=True)
sc.pl.umap(adata, color=['sample', 'total_counts', 'pct_counts_mt'])




# highly variable genes 선택
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, subset=True, batch_key="sampleID")

# scale
sc.pp.scale(adata, max_value=10)



n_batches = adata.var["highly_variable_nbatches"].value_counts()
ax = n_batches.plot(kind="bar")


# 1. 현재 상태의 adata (HVG 2000개만 포함)에서 복사
adata_scvi = adata.copy()

scvi.model.SCVI.setup_anndata(adata_scvi, layer="counts", batch_key="sampleID")
adata_scvi

import torch
print(torch.cuda.is_available())  # False이면 CPU로만 진행됨


model_scvi = scvi.model.SCVI(adata_scvi)
model_scvi

model_scvi.view_anndata_setup()

import numpy as np

max_epochs_scvi = np.min([round((20000 / adata.n_obs) * 400), 400])
max_epochs_scvi


model_scvi.train()

adata_scvi.obsm["X_scVI"] = model_scvi.get_latent_representation()

sc.pp.neighbors(adata_scvi, use_rep="X_scVI")
sc.tl.umap(adata_scvi)
adata_scvi

sc.pl.umap(adata_scvi,  color=['sample', 'pct_counts_mt', 'total_counts'], title='scVI UMAP')

from sklearn.metrics import silhouette_score

resolutions = [0.2, 0.4, 0.6, 0.8, 1.0]
scores = {}

for res in resolutions:
    sc.tl.leiden(adata_scvi, resolution=res, key_added=f'leiden_{res}')
    score = silhouette_score(adata_scvi.obsm['X_scVI'], adata_scvi.obs[f'leiden_{res}'])
    scores[res] = score
    print(f"🔹 resolution={res} → silhouette_score={score:.4f}")

# 최적 resolution 적용
best_res = max(scores, key=scores.get)
sc.tl.leiden(adata_scvi, resolution=best_res, key_added='leiden')
sc.pl.umap(adata_scvi, color='leiden', title=f"Leiden (res={best_res})")


# 🔹 클러스터별 marker gene 탐색
sc.tl.rank_genes_groups(
    adata_scvi,
    groupby='leiden_0.2',
    method='wilcoxon',  # 또는 't-test', 'logreg'
    key_added='rank_genes_leiden_0.6'
)

# 🔹 상위 marker 유전자 시각화
sc.pl.rank_genes_groups(
    adata_scvi,
    key='rank_genes_leiden_0.6',
    n_genes=5,  # 클러스터당 상위 5개
    sharey=False
)


# Celltypist용 AnnData 복사본 생성 (원본 유지)
adata_celltypist = adata_scvi.copy()

# raw count layer를 X로 설정
adata_celltypist.X = adata_celltypist.layers["counts"].copy()

# 세포 단위 정규화 (1만 counts로 맞춤)
import scanpy as sc
sc.pp.normalize_total(adata_celltypist, target_sum=1e4)

# log1p 변환
sc.pp.log1p(adata_celltypist)

# 희소행렬을 dense matrix로 변환 (celltypist는 dense 입력 요구)
import numpy as np
if not isinstance(adata_celltypist.X, np.ndarray):
    adata_celltypist.X = adata_celltypist.X.toarray()


from celltypist import models


models.download_models(
    force_update=True, 
    model=["Mouse_Whole_Brain.pkl"]
)

model_1 = models.Model.load(model="Mouse_Whole_Brain.pkl")

model_1.cell_types

import celltypist


predictions_high = celltypist.annotate(
    adata_celltypist, 
    model=model_1, 
    majority_voting=True
)

# 예측 라벨의 인덱스 중복 여부 확인
pred_idx = predictions_high.predicted_labels.index
duplicated = pred_idx.duplicated()
print(f"중복 인덱스 수: {duplicated.sum()}")
print(pred_idx[duplicated])


# 전체 셀 수
total_cells = len(predictions_high.predicted_labels)

# 중복 인덱스 개수
duplicated_cells = predictions_high.predicted_labels.index.duplicated(keep=False).sum()

# 중복 비율
dup_ratio = duplicated_cells / total_cells * 100

print(f"전체 셀 수: {total_cells}")
print(f"중복된 인덱스 수: {duplicated_cells}")
print(f"중복 비율: {dup_ratio:.2f}%")


# 예측 라벨의 인덱스 중복 여부 확인
pred_idx = predictions_high.predicted_labels.index
duplicated = pred_idx.duplicated()
print(f"중복 인덱스 수: {duplicated.sum()}")
print(pred_idx[duplicated])


# 중복된 index만 추출
duplicated_index = adata.obs_names[adata.obs_names.duplicated(keep=False)]

# 중복 인덱스에 해당하는 obs 데이터 확인
duplicated_obs = adata.obs.loc[duplicated_index]

# 중복 인덱스들의 sampleID를 집계
duplicated_obs['sampleID'].value_counts()


# 중복된 인덱스와 sampleID 매칭 확인
duplicated_obs[['sampleID']].groupby(duplicated_obs.index).value_counts()

# 중복된 index 중 첫 번째 것만 남기고 나머지는 제거
predictions_high.predicted_labels = predictions_high.predicted_labels[~predictions_high.predicted_labels.index.duplicated()]

predictions_high.to_adata()

# 예측 결과를 AnnData로 변환
predictions_high_adata = predictions_high.to_adata()

# index가 중복되지 않도록 정리 (예방 차원)
predictions_high_adata = predictions_high_adata[~predictions_high_adata.obs.index.duplicated(), :]

# celltypist 결과를 adata.obs에 추가
adata_scvi.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs.loc[adata.obs.index, "majority_voting"]
adata_scvi.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs.loc[adata.obs.index, "conf_score"]

adata_scvi.obs["celltypist_cell_label_coarse"] = predictions_high_adata.obs["majority_voting"]
adata_scvi.obs["celltypist_conf_score_coarse"] = predictions_high_adata.obs["conf_score"]



sc.pl.umap(
    adata_scvi,
    color=["celltypist_cell_label_coarse", "celltypist_conf_score_coarse"],
    frameon=False,
    sort_order=False,
    wspace=1,
)