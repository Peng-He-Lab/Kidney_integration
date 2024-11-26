import umap, scanpy as sc


for _path in ["Kidney_Subset_NoDuplicate_emb.h5ad", "Kidney_Subset_NoDuplicate_hvg5k.h5ad"]:

    adata = sc.read_h5ad(_path)
    obsm_list = list(adata.obsm.keys()).copy()

    for _obsm in obsm_list:
        if adata.obsm[_obsm].shape[1] > 2:
            print(_obsm)
            umap_transformer = umap.UMAP(min_dist=0.5)
            adata.obsm[_obsm + "_umap"] = umap_transformer.fit_transform(adata.obsm[_obsm])

    adata.write(_path, compression="gzip")
