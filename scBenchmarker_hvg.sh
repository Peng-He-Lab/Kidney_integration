# python scBenchmarker.py \
#     --adata_path Kidney_Combined_hvg5k.h5ad \
#     --save_path Kidney_Combined_hvg5k.h5ad \
#     --highvar --highvar_n 5000 \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --gpu 0 \
#     --savecsv kidney_hvg5k \
#     --saveadata --all;

# Kidney_Combined_cellbender_v2.h5ad

python scBenchmarker.py \
    --adata_path Kidney_Combined_hvg5k_LabeledSubset.h5ad \
    --save_path Kidney_Combined_hvg5k_LabeledSubset.h5ad \
    --batch_key batch \
    --gpu 0 \
    --label_key harmonized_celltype \
    --savecsv kidney_hvg5k \
    --saveadata --all;
