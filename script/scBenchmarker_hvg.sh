cd ..

# python src/scBenchmarker.py \
#     --adata_path data/Kidney_Combined_cellbender_v2.h5ad \
#     --save_path data/Kidney_Combined_hvg15k.h5ad \
#     --highvar \
#     --highvar_n 15885 \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --savecsv kidney_hvg15k \
#     --saveadata \
#     --all;


python src/labeled_subset.py \
    data/Kidney_Combined_hvg15k.h5ad \
    data/Kidney_Combined_hvg15k_LabeledSubset.h5ad

python src/scBenchmarker.py \
    --adata_path data/Kidney_Combined_hvg15k_LabeledSubset.h5ad \
    --save_path data/Kidney_Combined_hvg15k_LabeledSubset.h5ad \
    --batch_key batch \
    --label_key harmonized_celltype \
    --savecsv kidney_hvg15k \
    --saveadata \
    --all;


# Kidney_Combined_cellbender_v2.h5ad

# python scBenchmarker.py \
#     --adata_path Kidney_Combined_hvg5k_LabeledSubset.h5ad \
#     --save_path Kidney_Combined_hvg5k_LabeledSubset.h5ad \
#     --batch_key batch \
#     --gpu 0 \
#     --label_key harmonized_celltype \
#     --savecsv kidney_hvg5k \
#     --saveadata --all;
