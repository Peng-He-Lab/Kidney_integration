cd ..

# python src/scBenchmarker.py \
#     --adata_path data/Kidney_Combined_deeptree.h5ad \
#     --save_path data/Kidney_Combined_deeptree.h5ad \
#     --batch_key batch \
#     --label_key harmonized_celltype \
#     --savecsv kidney_deeptree \
#     --saveadata \
#     --all;

python src/labeled_subset.py \
    data/Kidney_Combined_deeptree.h5ad \
    data/Kidney_Combined_deeptree_LabeledSubset.h5ad

python src/scBenchmarker.py \
    --adata_path data/Kidney_Combined_deeptree_LabeledSubset.h5ad \
    --save_path data/Kidney_Combined_deeptree_LabeledSubset.h5ad \
    --batch_key batch \
    --label_key harmonized_celltype \
    --savecsv kidney_deeptree \
    --saveadata \
    --all;

# python src/umap.py