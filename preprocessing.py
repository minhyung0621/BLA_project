import scanpy as sc
import os
from pathlib import Path

# ✅ 1. 경로 설정 및 초기 폴더 생성
base_data_dir = "/data1/BLA_projection_2/Data/Resources/.v2/filtered-feature_bc_matrix"  # 여러 샘플이 있는 디렉토리
save_h5ad_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_adata생성"        # h5ad 저장 경로
save_fig_dir = "/data1/BLA_projection_2/Data/Resources/.v2/figure"              # 그림 저장 경로

os.makedirs(save_h5ad_dir, exist_ok=True)
os.makedirs(save_fig_dir, exist_ok=True)

# ✅ 2. 모든 샘플 디렉토리 순회
for sample_dir in Path(base_data_dir).iterdir():
    if not sample_dir.is_dir():
        continue
    mtx_path = sample_dir / "filtered_feature_bc_matrix"
    if not mtx_path.exists():
        continue  # mtX 디렉토리 없으면 건너뜀

    sample_name = f"BLA_projection_sample_{sample_dir.name}"  # 샘플명 지정

    # ✅ 3. adata 로딩 및 QC 메트릭 계산
    adata = sc.read_10x_mtx(mtx_path, var_names='gene_symbols', cache=True)
    adata.obs["sampleID"] = sample_name
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # ✅ 4. .h5ad 파일 저장
    h5ad_save_path = os.path.join(save_h5ad_dir, f"{sample_name}_after_adata생성.h5ad")
    adata.write(h5ad_save_path)  # # 각 샘플별 h5ad 파일 저장

    # ✅ 5. QC plot 저장
    sc.settings.figdir = save_fig_dir  # # 그림 저장 폴더 지정
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        color="sampleID",
        title=f"{sample_name}: Total counts vs Mitochondrial %",
        save=f"{sample_name}_total_vs_mito_pct.png"  # # 각 샘플별로 그림 이름 다르게 저장
    )

print("✅ 모든 샘플에 대한 QC 및 저장이 완료되었습니다.")


import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import os

# 🔹 경로 설정
input_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_adata생성"
violin_plot_dir = "/data1/BLA_projection_2/Data/Resources/.v2/figure/after_QC_violin_plot"
after_qc_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_QC"

os.makedirs(violin_plot_dir, exist_ok=True)
os.makedirs(after_qc_dir, exist_ok=True)

# 🔹 QC 지표 설정
qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]

# 🔹 adata 파일 목록
adata_files = [f for f in os.listdir(input_dir) if f.endswith(".h5ad")]
adata_list = {}

for file in adata_files:
    sample_name = os.path.splitext(file)[0]
    file_path = os.path.join(input_dir, file)
    adata = sc.read_h5ad(file_path)
    adata_list[sample_name] = adata

# 🔹 전체 violin plot 그리기
n_metrics = len(qc_metrics)
n_samples = len(adata_list)
fig, axs = plt.subplots(
    nrows=n_samples,
    ncols=n_metrics,
    figsize=(4 * n_metrics, 3 * n_samples),
    squeeze=False
)

summary_stats = []

for row_idx, (sample, adata) in enumerate(adata_list.items()):
    sample_id = list(adata.obs['sampleID'].unique())[0] if 'sampleID' in adata.obs.columns else sample
    for col_idx, metric in enumerate(qc_metrics):
        ax = axs[row_idx, col_idx]
        sc.pl.violin(adata, keys=metric, jitter=0.4, ax=ax, show=False)
        if col_idx == 0:
            ax.set_title(sample_id, fontsize=12, loc='left')
        else:
            ax.set_ylabel('')
        values = adata.obs[metric]
        stats = values.describe()
        stats['median'] = values.median()
        stats_dict = stats.to_dict()
        stats_dict['sample'] = sample
        stats_dict['metric'] = metric
        summary_stats.append(stats_dict)

# 🔹 전체 sample 묶은 plot 저장
plt.tight_layout()
combined_plot_path = os.path.join(violin_plot_dir, "ALL_SAMPLES_QC_violin_combined.pdf")
plt.savefig(combined_plot_path, bbox_inches='tight')
plt.show()

# 🔹 개별 violin plot & QC 후 파일 저장
for sample, adata in adata_list.items():
    # violin plot 저장
    fig, axs = plt.subplots(1, len(qc_metrics), figsize=(4 * len(qc_metrics), 4))
    for i, metric in enumerate(qc_metrics):
        sc.pl.violin(adata, keys=metric, jitter=0.4, ax=axs[i], show=False)
        axs[i].set_title(metric)
    plt.tight_layout()
    fig.savefig(os.path.join(violin_plot_dir, f"{sample}_QC_violin.pdf"))
    plt.close(fig)

    # QC 후 파일 저장
    adata.write(os.path.join(after_qc_dir, f"{sample}_after_QC.h5ad"))

# 🔹 통계 요약 화면 출력
qc_summary_df = pd.DataFrame(summary_stats)
qc_summary_df = qc_summary_df[[ 'sample', 'metric', 'mean', 'std', 'min', '25%', 'median', '50%', '75%', 'max' ]]
print(qc_summary_df)


# 🔹 통계 요약 가공 (sample별, metric별로 pivot)
qc_summary_df = pd.DataFrame(summary_stats)

# 필요한 통계만 선택
selected_stats = ['mean', 'median']
wide_format_rows = []

for sample in qc_summary_df['sample'].unique():
    sample_df = qc_summary_df[qc_summary_df['sample'] == sample]
    row = {'sample': sample}
    for _, row_data in sample_df.iterrows():
        metric = row_data['metric']
        for stat in selected_stats:
            row[f'{metric}_{stat}'] = row_data[stat]
    wide_format_rows.append(row)

# DataFrame 생성 및 정렬
final_df = pd.DataFrame(wide_format_rows)
final_df = final_df.sort_values(by='sample')

# 출력
pd.set_option('display.max_columns', None)  # 모든 컬럼 표시
print(final_df)



import os
import pandas as pd
import scrublet as scr
import scanpy as sc
import anndata as ad

def run_preprocess(adata_tmp):
    adata_tmp = adata_tmp.copy()
    sample_id = list(adata_tmp.obs['sampleID'].unique())[0]
    print(f'Apply QC for {sample_id}')

    # Ensure var_names are proper gene identifiers
    if not pd.api.types.is_string_dtype(adata_tmp.var_names):
        if 'gene_symbols' in adata_tmp.var.columns:
            adata_tmp.var_names = adata_tmp.var['gene_symbols']
        elif 'gene_ids' in adata_tmp.var.columns:
            adata_tmp.var_names = adata_tmp.var['gene_ids']
        adata_tmp.var_names_make_unique()
        adata_tmp.var.index = adata_tmp.var_names

    # Basic filtering
    sc.pp.filter_genes(adata_tmp, min_cells=3)
    adata_tmp = adata_tmp[
        (adata_tmp.obs.n_genes_by_counts >= 500) &
        (adata_tmp.obs.n_genes_by_counts <= 7000) &
        (adata_tmp.obs.total_counts >= 1000) &
        (adata_tmp.obs.pct_counts_mt < 20),
        :
    ]

    # Scrublet filtering
    try:
        scrub = scr.Scrublet(adata_tmp.X)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False, n_prin_comps=20)
        adata_tmp.obs['doublet_scores'] = doublet_scores
        adata_tmp.obs['predicted_doublets'] = predicted_doublets.astype(str)
        adata_tmp = adata_tmp[adata_tmp.obs.doublet_scores < 0.2, :]
    except Exception as e:
        print(f'#Pass {sample_id}: Scrublet failed ({e})')
        return None

    if adata_tmp.n_obs < 500:
        print(f"#Skip {sample_id}: Less than 500 cells after QC")
        return None

    print(f'✅ Preprocessing complete for {sample_id}')
    return adata_tmp

# ▶ 분석할 h5ad 파일들이 위치한 디렉토리
input_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_QC"
output_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_preprocessing"
os.makedirs(output_dir, exist_ok=True)

# ▶ 파일 목록 가져오기
file_list = [f for f in os.listdir(input_dir) if f.endswith(".h5ad")]
file_list.sort()  # 정렬 (선택사항)

# ▶ 각 파일에 대해 QC + 저장 수행
for filename in file_list:
    input_path = os.path.join(input_dir, filename)
    print(f"\n📂 Processing: {filename}")
    
    try:
        adata = sc.read_h5ad(input_path)
        adata_qc = run_preprocess(adata)

        if adata_qc is not None:
            sample_id = list(adata_qc.obs['sampleID'].unique())[0]
            out_file = f"{sample_id}_after_preprocessing.h5ad"
            out_path = os.path.join(output_dir, out_file)
            adata_qc.write(out_path)
            print(f"✅ Saved to {out_path}")
        else:
            print("❌ QC 실패: 저장하지 않음 (None 반환됨)")
    except Exception as e:
        print(f"🚫 Failed to process {filename}: {e}")


import matplotlib.pyplot as plt
from PIL import Image
import os

# 입력 디렉토리 경로
input_dir = "/data1/BLA_projection_2/Data/Resources/.v2/figure/after_adata_generation_mt"

# 파일 목록 가져오기 (확장자가 .png인 것만)
image_files = sorted([f for f in os.listdir(input_dir) if f.endswith(".png")])

# 이미지 불러오기
images = [Image.open(os.path.join(input_dir, f)) for f in image_files]

# 전체 figure 설정 (2행 3열로 배치)
fig, axs = plt.subplots(2, 3, figsize=(18, 10))
axs = axs.flatten()

# 각 subplot에 이미지 삽입
for ax, img, name in zip(axs, images, image_files):
    ax.imshow(img)
    ax.axis('off')
    ax.set_title(name.replace(".png", ""), fontsize=10)

# 전체 제목 및 저장
plt.suptitle("Total Counts vs Mitochondrial % Across Samples", fontsize=16)
plt.tight_layout(rect=[0, 0, 1, 0.95])

# 저장 경로
out_path = os.path.join(input_dir, "merged_total_vs_mt.png")
plt.savefig(out_path, dpi=300)
plt.show()

print(f"✅ 병합된 이미지 저장 완료: {out_path}")



import os
import scanpy as sc
import pandas as pd

# 경로 지정
data_dir = "/data1/BLA_projection_2/Data/Resources/.v2/after_preprocessing"

# h5ad 파일 목록 가져오기
h5ad_files = [f for f in os.listdir(data_dir) if f.endswith(".h5ad")]

# 결과 저장용 리스트
summary_list = []

# 각 파일에 대해 gene/cell 수 계산
for file in h5ad_files:
    file_path = os.path.join(data_dir, file)
    adata = sc.read_h5ad(file_path)
    sample_id = file.replace(".h5ad", "")
    
    cell_count = adata.n_obs
    gene_count = adata.n_vars

    summary_list.append({
        "Sample": sample_id,
        "Cell Count": cell_count,
        "Gene Count": gene_count
    })

# 결과를 DataFrame으로 출력
summary_df = pd.DataFrame(summary_list)
print(summary_df)

