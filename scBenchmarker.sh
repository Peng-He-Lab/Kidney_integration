CUDA_VISIBLE_DEVICES=7 & \
python scBenchmarker.py \
    --adata_path Kidney_Subset_DeepTree15k.h5ad \
    --save_path Kidney_Subset_DeepTree15k.h5ad \
    --batch_key batch \
    --label_key harmonized_celltype \
    --savecsv kidney_deeptree \
    --saveadata --all;

# python scBenchmarker.py \
#     --adata_path Kidney_Subset_emb.h5ad \
#     --save_path Kidney_Subset_emb.h5ad \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --all;
