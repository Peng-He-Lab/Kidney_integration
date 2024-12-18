# python scBenchmarker.py \
#     --adata_path Kidney_Combined_sn_hvg5k_LabeledSubset.h5ad \
#     --save_path Kidney_Combined_sn_hvg5k_LabeledSubset.h5ad \
#     --highvar --highvar_n 5000 \
#     --gpu 2 \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --savecsv kidney_hvg5k_sn \
#     --saveadata --all;

cd ..

python src/scBenchmarker.py \
    --adata_path data/Kidney_Combined_sc_hvg5k_LabeledSubset.h5ad \
    --save_path data/Kidney_Combined_sc_hvg5k_LabeledSubset.h5ad \
    --highvar \
    --highvar_n 5000 \
    --batch_key batch \
    --gpu 3 \
    --label_key harmonized_celltype \
    --savecsv kidney_hvg5k_sc \
    --saveadata \
    --all;


# python scBenchmarker.py \
#     --adata_path Kidney_Combined_sn_hvg5k.h5ad \
#     --save_path Kidney_Combined_sn_hvg5k.h5ad \
#     --highvar --highvar_n 5000 \
#     --gpu 2 \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --savecsv kidney_hvg5k_sn \
#     --saveadata --all;

# python scBenchmarker.py \
#     --adata_path Kidney_Combined_sc_hvg5k.h5ad \
#     --save_path Kidney_Combined_sc_hvg5k.h5ad \
#     --highvar --highvar_n 5000 \
#     --batch_key batch \
#     --gpu 3 \
#     --label_key harmonized_celltype \
#     --savecsv kidney_hvg5k_sc \
#     --saveadata --all;
