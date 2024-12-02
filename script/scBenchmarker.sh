cd ..

# CUDA_VISIBLE_DEVICES=7 & \
# python src/scBenchmarker.py \
#     --adata_path data/Kidney_Combined_whole.h5ad  \
#     --save_path data/Kidney_Combined_whole.h5ad \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --savecsv kidney_whole \
#     --saveadata \
#     --all;

python src/labeled_subset.py \
    data/Kidney_Combined_whole.h5ad \
    data/Kidney_Combined_whole_LabeledSubset.h5ad

python src/scBenchmarker.py \
    --adata_path data/Kidney_Combined_whole_LabeledSubset.h5ad \
    --save_path data/Kidney_Combined_whole_LabeledSubset.h5ad \
    --batch_key batch \
    --label_key harmonized_celltype \
    --savecsv kidney_whole \
    --saveadata \
    --all;


# python scBenchmarker.py \
#     --adata_path Kidney_Subset_emb.h5ad \
#     --save_path Kidney_Subset_emb.h5ad \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --all;
