import umap, scanpy as sc

# for _file in ["Kidney_Subset_NoDuplicate_emb.h5ad", "Kidney_Subset_NoDuplicate_hvg5k.h5ad"]:
for _file in [
    "Kidney_Combined_sn_hvg5k.h5ad", 
    "Kidney_Combined_sn_hvg5k_LabeledSubset.h5ad",
    "Kidney_Combined_sc_hvg5k.h5ad", 
    "Kidney_Combined_sc_hvg5k_LabeledSubset.h5ad",
    "Kidney_Combined_hvg5k.h5ad", 
    "Kidney_Combined_hvg5k_LabeledSubset.h5ad"
    ]:
    print(_file)
    adata = sc.read_h5ad(_file)
    obsm_list = list(adata.obsm.keys()).copy()

    for _obsm in obsm_list:
        if adata.obsm[_obsm].shape[1] > 2:
            if _obsm + "_umap" in adata.obsm.keys():
                continue
            print(_obsm)
            umap_transformer = umap.UMAP(min_dist=0.5)
            adata.obsm[_obsm + "_umap"] = umap_transformer.fit_transform(adata.obsm[_obsm])

    adata.write(_file, compression="gzip")
    print("=" * 77)
