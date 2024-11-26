# https://docs.scvi-tools.org/en/stable/tutorials/notebooks/harmonization.html
# https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html

import warnings
warnings.filterwarnings("ignore")

import os, scvi, torch, argparse, numpy as np, scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
# from os.path import join, exists, getmtime, dirname
# from torch.utils.data import DataLoader


# set_seed(dh.SEED)
# os.environ["CUDA_VISIBLE_DEVICES"] = dh.CUDA_DEVICE
# torch.cuda.is_available() -> False
sc.settings.verbosity = 3
sc.logging.print_header()

GPU_4_NEIGHBORS = False

print("\n\n")
print("=" * 44)
print("use GPU for neighbors calculation: ", GPU_4_NEIGHBORS)
print("GPU is available: ", torch.cuda.is_available())
print("=" * 44)
print("\n")

RES_DIR = "/scratch/site/u/wangh256/kidney_atlas/res"

def faiss_hnsw_nn(X: np.ndarray, k: int):
    # """GPU HNSW nearest neighbor search using faiss.

    # See https://github.com/nmslib/hnswlib/blob/master/ALGO_PARAMS.md
    # for index param details.
    # """
    # X = np.ascontiguousarray(X, dtype=np.float32)
    # res = faiss.StandardGpuResources()
    # M = 32
    # index = faiss.IndexHNSWFlat(X.shape[1], M, faiss.METRIC_L2)
    # gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    # gpu_index.add(X)
    # distances, indices = gpu_index.search(X, k)
    # del index
    # del gpu_index
    # # distances are squared
    # return NeighborsOutput(indices=indices, distances=np.sqrt(distances))
    raise NotImplementedError


def faiss_brute_force_nn(X: np.ndarray, k: int):
    # """GPU brute force nearest neighbor search using faiss."""
    # X = np.ascontiguousarray(X, dtype=np.float32)
    # res = faiss.StandardGpuResources()
    # index = faiss.IndexFlatL2(X.shape[1])
    # gpu_index = faiss.index_cpu_to_gpu(res, 0, index)
    # gpu_index.add(X)
    # distances, indices = gpu_index.search(X, k)
    # del index
    # del gpu_index
    # # distances are squared
    # return NeighborsOutput(indices=indices, distances=np.sqrt(distances))
    raise NotImplementedError


def Parser_Benchmarker():
    parser = argparse.ArgumentParser(description="Parser for scBenchmarker")

    parser.add_argument("--gpu", type=str, default="0")
    # parser.add_argument("--save_path", type=str, default=rf"{dh.MODEL_DIR}/_minimal_")
    # parser.add_argument("--batch_id", type=str, nargs="+", default=None)
    # parser.add_argument("--dataset", type=str)
    parser.add_argument("--highvar", action="store_true")
    parser.add_argument("--highvar_n", type=int, default=5000)
    parser.add_argument("--savecsv", type=str, default=None)

    parser.add_argument("--adata_path", type=str)
    parser.add_argument("--save_path", type=str)
    parser.add_argument("--batch_key", type=str)
    parser.add_argument("--label_key", type=str)

    # === Integration Baselines ===
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--umap", action="store_true")
    parser.add_argument("--tsne", action="store_true")
    parser.add_argument("--scvi", action="store_true")
    parser.add_argument("--bbknn", action="store_true")
    parser.add_argument("--scgen", action="store_true")
    parser.add_argument("--scpoli", action="store_true")
    parser.add_argument("--scanvi", action="store_true")
    parser.add_argument("--harmony", action="store_true")
    parser.add_argument("--fastmnn", action="store_true")
    parser.add_argument("--scanorama", action="store_true")
    #
    # parser.add_argument("--islander", action="store_true")
    parser.add_argument("--obsm_keys", type=str, nargs="+", default=None)

    parser.add_argument("--use_raw", action="store_true")
    parser.add_argument("--saveadata", action="store_true")

    #  === local debug ===
    parser.add_argument("--debug", action="store_true", default=False)
    return parser.parse_args()



class scIB():
    """Integration Benchmarker"""

    # https://scib-metrics.readthedocs.io/en/stable/
    # https://github.com/theislab/scib-reproducibility/
    # https://github.com/theislab/scib-pipeline/
    # https://github.com/theislab/scib/

    def __init__(self, args):
        self.args = args
        # self.device = torch.device("cuda:%s" % self.args.gpu if torch.cuda.is_available() else "cpu")
        # os.environ["CUDA_VISIBLE_DEVICES"] = self.args.gpu
        self._load_adata_()

    def _load_adata_(self):
        
        assert not self.args.use_raw, "use_raw is not supported"
        self.adata = sc.read(self.args.adata_path)
        self.batch_key = self.args.batch_key
        self.label_key = self.args.label_key
        
        # check if adata.X is raw counts
        sc.pp.filter_cells(self.adata, min_genes=1)
        sc.pp.filter_genes(self.adata, min_counts=1)

        subset = self.adata.X[:100, :].toarray()
        non_negative = np.all(subset >= 0)
        integer_values = np.all(subset.astype(int) == subset)
        assert non_negative  
        
        if not integer_values:
            self._str_formatter("reset the counts to adata.X, then norm + log1p")
            self.adata.X = self.adata.layers["counts"].copy()
        else:
            self._str_formatter("reset the adata.X to counts, then norm + log1p")
            self.adata.layers["counts"] = self.adata.X.copy()
        
        sc.pp.normalize_total(self.adata, target_sum=1e4)
        sc.pp.log1p(self.adata)
        self.adata.layers["log1p_norm"] = self.adata.X

        if self.args.highvar:
            try:
                sc.pp.highly_variable_genes(
                    self.adata, 
                    flavor="seurat_v3", 
                    n_top_genes=self.args.highvar_n, 
                    batch_key=self.batch_key,
                    subset=True
                    )
                # sc.pp.highly_variable_genes(
                #     self.adata, 
                #     flavor="seurat_v3_paper", 
                #     n_top_genes=self.args.highvar_n, 
                #     batch_key=self.batch_key
                #     )
                message = f"select {self.args.highvar_n} hvg genes, batchwise: {self.batch_key}"

            except:
                sc.pp.highly_variable_genes(
                    self.adata, 
                    flavor="seurat_v3", 
                    n_top_genes=self.args.highvar_n, 
                    subset=True
                    )
                message = f"select {self.args.highvar_n} hvg genes, overall"
            
            self._str_formatter(message)
            
        return

    def _str_formatter(self, message):
        print(f"\n=== {message} ===\n")

    def _scDataloader(self):
        # self._str_formatter("Dataloader")
        # _verb = True
        # _scDataset = scDataset(
        #     dataset=self.args.dataset,
        #     inference=True,
        #     rm_cache=False,
        #     verbose=_verb,
        # )
        # _scDataLoader = DataLoader(
        #     _scDataset,
        #     batch_size=1,
        #     shuffle=False,
        #     num_workers=8,
        #     collate_fn=collate_fn,
        # )
        # return _scDataset, _scDataLoader
        return None

    def _benchmark_(self, n_jobs=-1, scratch=False):
        # NOTE: for DEBUG
        # first_5_ = self.adata.obs[self.batch_key].unique()[:5].to_list()
        # self.adata = self.adata[[_it in first_5_ for _it in self.adata.obs[self.batch_key]]]

        # recompute the embeddings
        # self.scratch = True if self.args.highvar else scratch
        self.scratch = scratch
        self._pca_()
        if self.args.umap or self.args.all:
            self._umap_()
        if self.args.tsne or self.args.all:
            self._tsne_()
        if self.args.harmony or self.args.all:
            self._harmony_()
        if self.args.scanorama or self.args.all:
            self._scanorama_()
        if self.args.scvi or self.args.all:
            self._scvi_()
        if self.args.scanvi or self.args.all:
            self._scanvi_()
        if self.args.bbknn or self.args.all:
            self._bbknn_()
        if self.args.scgen or self.args.all:
            self._scgen_()
        if self.args.fastmnn or self.args.all:
            self._fastmnn_()
        if self.args.scpoli or self.args.all:
            self._scpoli_()

        # if self.args.islander:
        #     self._islander_()

        if self.args.obsm_keys is None:
            obsm_keys = list(self.adata.obsm)
        else:
            obsm_keys = self.args.obsm_keys

        for embed in obsm_keys:
            print("%12s, %d" % (embed, self.adata.obsm[embed].shape[1]))
            if self.adata.obsm[embed].shape[0] != np.unique(self.adata.obsm[embed], axis=0).shape[0]:
                # print("\nWarning: Embedding %s has duplications, Dropped.\n" % embed)
                print("\nWarning: Embedding %s has duplications.\n" % embed)
                # obsm_keys.remove(embed)
        # self._save_adata_()

        self._str_formatter(rf"scIB Benchmarking: {obsm_keys}")
        biocons = BioConservation(nmi_ari_cluster_labels_leiden=True)
        batcorr = BatchCorrection()

        """ === for DEBUG ==="""
        # biocons = BioConservation(isolated_labels=False)
        # biocons = BioConservation(
        #     silhouette_label=False,
        #     isolated_labels=False,)
        # batcorr = BatchCorrection(
        #     silhouette_batch=False,
        #     ilisi_knn=False,
        #     kbet_per_label=False,
        # )

        self.benchmarker = Benchmarker(
            self.adata,
            batch_key=self.batch_key,
            label_key=self.label_key,
            embedding_obsm_keys=obsm_keys,
            bio_conservation_metrics=biocons,
            batch_correction_metrics=batcorr,
            pre_integrated_embedding_obsm_key="X_pca",
            n_jobs=n_jobs,
        )
        if torch.cuda.is_available() and GPU_4_NEIGHBORS:
            # self.benchmarker.prepare(neighbor_computer=faiss_brute_force_nn)
            self.benchmarker.prepare(neighbor_computer=faiss_hnsw_nn)
        else:
            # Calculate the Neighbors based on the CPUs
            self.benchmarker.prepare(neighbor_computer=None)
        self.benchmarker.benchmark()

        self._str_formatter("scIB Benchmarking Finished")
        df = self.benchmarker.get_results(min_max_scale=False)
        print(df.head(7))

        os.makedirs(rf"./scIB", exist_ok=True)
        if self.args.savecsv is not None:
            df.to_csv(rf"./scIB/{self.args.savecsv}.csv")
        else:
            df.to_csv(rf"./scIB/result.csv")

        # savefig = False
        # if savefig:
        #     _suffix = "_hvg" if self.args.highvar else ""
        #     self.benchmarker.plot_results_table(min_max_scale=False, show=False, save_dir=rf"{dh.RES_DIR}/scIB/")
        #     os.rename(
        #         src=rf"{dh.RES_DIR}/scIB/scib_results.svg",
        #         dst=rf"{dh.RES_DIR}/figures/scib_{self.args.dataset}{_suffix}.svg",
        #     )
        return

    def _save_adata_(self):
        if self.args.saveadata:
            if self.args.save_path is not None:
                _save_path = self.args.save_path
            else:
                _save_path = self.args.adata_path.replace(".h5ad", "_emb.h5ad")
            self._str_formatter("Saving %s" % _save_path)
            self.adata.write_h5ad(_save_path, compression="gzip")
        return

    def _pca_(self, n_comps=50):
        self._str_formatter("PCA")
        if "X_pca" in self.adata.obsm and not self.scratch:
            return
        try:
            sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3", n_top_genes=self.args.highvar_n, batch_key=self.batch_key)
        except:
            sc.pp.highly_variable_genes(self.adata, flavor="seurat_v3", n_top_genes=self.args.highvar_n)
        sc.pp.pca(self.adata, n_comps=n_comps, use_highly_variable=True, svd_solver="arpack")
        self._save_adata_()
        return

    def _umap_(self, n_neighbors=10, n_pcs=50):
        self._str_formatter("UMAP")
        if "X_umap" in self.adata.obsm and not self.scratch:
            return
        sc.pp.neighbors(self.adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(self.adata)
        self._save_adata_()
        return

    def _tsne_(self):
        if "X_tsne" in self.adata.obsm and not self.scratch:
            return
        self._str_formatter("TSNE")
        sc.tl.tsne(self.adata)
        self._save_adata_()
        return

    def _scvi_(self):
        self._str_formatter("scVI")
        if "X_scVI" in self.adata.obsm and not self.scratch:
            return
        adata = self.adata.copy()
        adata.X = adata.layers["counts"]
        
        scvi.model.SCVI.setup_anndata(adata, layer=None, batch_key=self.batch_key)
        self.vae = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
        self.vae.train()
        self.adata.obsm["X_scVI"] = self.vae.get_latent_representation()
        self._save_adata_()
        return

    def _scanvi_(self):
        self._str_formatter("X_scANVI")
        if "X_scANVI" in self.adata.obsm and not self.scratch:
            return
        try:
            lvae = scvi.model.SCANVI.from_scvi_model(
                self.vae,
                adata=self.adata,
                labels_key=self.label_key,
                unlabeled_category="Unknown",
            )
            lvae.train()
            self.adata.obsm["X_scANVI"] = lvae.get_latent_representation(self.adata)
            self._save_adata_()
        except:
            self._str_formatter("scANVI skipped")
        return

    def _bbknn_(self):
        # ref: https://bbknn.readthedocs.io/en/latest/#
        # tutorial: https://nbviewer.org/github/Teichlab/bbknn/blob/master/examples/mouse.ipynb

        import bbknn

        self._str_formatter("bbknn")
        if "X_bbknn" in self.adata.obsm and not self.scratch:
            return
        # if self.adata.n_obs < 1e5:
        #     _temp_adata = bbknn.bbknn(self.adata, batch_key=self.batch_key, copy=True)
        # else:
        print(self.adata.obs[self.batch_key].value_counts().tail())
        _smallest_n_neighbor = self.adata.obs[self.batch_key].value_counts().tail(1).values[0]
        _temp_adata = bbknn.bbknn(
            self.adata,
            batch_key=self.batch_key,
            neighbors_within_batch=min(10, _smallest_n_neighbor),
            copy=True,
        )
        sc.tl.umap(_temp_adata)
        self.adata.obsm["X_bbknn"] = _temp_adata.obsm["X_umap"]
        self._save_adata_()
        return

    def _harmony_(self):
        # https://pypi.org/project/harmony-pytorch/
        self._str_formatter("Harmony")
        if "Harmony" in self.adata.obsm and not self.scratch:
            return

        from harmony import harmonize

        self.adata.obsm["Harmony"] = harmonize(
            self.adata.obsm["X_pca"], 
            self.adata.obs, 
            batch_key=self.batch_key, 
            use_gpu=True
            )
        self._save_adata_()
        return

    def _scanorama_(self):
        # https://github.com/brianhie/scanorama
        self._str_formatter("Scanorama")
        if "Scanorama" in self.adata.obsm and not self.scratch:
            return

        import scanorama

        batch_cats = self.adata.obs[self.batch_key].cat.categories
        adata_list = [self.adata[self.adata.obs[self.batch_key] == b].copy() for b in batch_cats]
        scanorama.integrate_scanpy(adata_list)

        self.adata.obsm["Scanorama"] = np.zeros((self.adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))
        for i, b in enumerate(batch_cats):
            self.adata.obsm["Scanorama"][self.adata.obs[self.batch_key] == b] = adata_list[i].obsm["X_scanorama"]
        self._save_adata_()
        return

    def _scgen_(self):
        self._str_formatter("scGen")
        if "scgen_pca" in self.adata.obsm and not self.scratch:
            return

        # if not self.args.highvar:
        if self.adata.shape[0] > 1e5:
            return

        import scgen

        # ref: https://scgen.readthedocs.io/en/stable/tutorials/scgen_batch_removal.html
        # pip install git+https://github.com/theislab/scgen.git

        # ref: https://github.com/theislab/scib-reproducibility/blob/main/notebooks/integration/Test_scgen.ipynb
        # ref: https://github.com/theislab/scib/blob/main/scib/integration.py

        # ref: https://github.com/LungCellAtlas/HLCA_reproducibility/blob/main/notebooks/1_building_and_annotating_the_atlas_core/04_integration_benchmark_prep_and_scgen.ipynb

        scgen.SCGEN.setup_anndata(self.adata, batch_key=self.batch_key, labels_key=self.label_key)
        model = scgen.SCGEN(self.adata)
        model.train(
            max_epochs=100,
            batch_size=32,  # 32
            early_stopping=True,
            early_stopping_patience=25,
        )

        adata_scgen = model.batch_removal()
        sc.pp.pca(adata_scgen, svd_solver="arpack")
        sc.pp.neighbors(adata_scgen)
        sc.tl.umap(adata_scgen)

        self.adata.obsm["scgen_umap"] = adata_scgen.obsm["X_umap"]
        self.adata.obsm["scgen_pca"] = adata_scgen.obsm["X_pca"]
        self._save_adata_()
        return

    def _fastmnn_(self):
        # [deprecated]: https://github.com/chriscainx/mnnpy
        # ref: https://github.com/HelloWorldLTY/mnnpy
        # from: https://github.com/chriscainx/mnnpy/issues/42

        self._str_formatter("fastMNN")
        if "fastMNN_pca" in self.adata.obsm and not self.scratch:
            return

        # if not self.args.highvar:
        if self.adata.shape[0] > 1e5:
            return

        import mnnpy

        def split_batches(adata, batch_key, hvg=None, return_categories=False):
            """Split batches and preserve category information
            Ref: https://github.com/theislab/scib/blob/main/scib/utils.py#L32"""
            split = []
            batch_categories = adata.obs[batch_key].cat.categories
            if hvg is not None:
                adata = adata[:, hvg]
            for i in batch_categories:
                split.append(adata[adata.obs[batch_key] == i].copy())
            if return_categories:
                return split, batch_categories
            return split

        split, categories = split_batches(self.adata, batch_key=self.batch_key, return_categories=True)
        # if self.args.dataset in [
        #     "lung_fetal_organoid",
        #     "COVID",
        #     "heart",
        #     "brain",
        #     "breast",
        # ]:
        #     k = 10
        # else:
        #     k = 20
        k = 20
        corrected, _, _ = mnnpy.mnn_correct(
            *split,
            k=k,
            batch_key=self.batch_key,
            batch_categories=categories,
            index_unique=None,
        )

        adata_fastmnn = corrected
        sc.pp.pca(adata_fastmnn, svd_solver="arpack")
        sc.pp.neighbors(adata_fastmnn)
        sc.tl.umap(adata_fastmnn)

        self.adata.obsm["fastMNN_umap"] = adata_fastmnn.obsm["X_umap"]
        self.adata.obsm["fastMNN_pca"] = adata_fastmnn.obsm["X_pca"]
        self._save_adata_()
        return

    def _scpoli_(self):
        self._str_formatter("scPoli")
        if "scPoli" in self.adata.obsm and not self.scratch:
            return
        # if self.args.dataset in ["brain", "breast"]:
        #     return
        from scarches.models.scpoli import scPoli

        self.adata.X = self.adata.X.astype(np.float32)
        early_stopping_kwargs = {
            "early_stopping_metric": "val_prototype_loss",
            "mode": "min",
            "threshold": 0,
            "patience": 10,
            "reduce_lr": True,
            "lr_patience": 13,
            "lr_factor": 0.1,
        }
        scpoli_model = scPoli(
            adata=self.adata,
            condition_keys=self.batch_key,
            cell_type_keys=self.label_key,
            embedding_dims=16,
            recon_loss="nb",
        )
        scpoli_model.train(
            n_epochs=50,
            pretraining_epochs=40,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=5,
        )
        self.adata.obsm["scPoli"] = scpoli_model.get_latent(self.adata, mean=True)
        self._save_adata_()
        return

    def _islander_(self):
        # from tqdm.auto import tqdm
        # import umap
        # scDataset, self.scDataLoader = self._scDataloader()
        # # self.cell2cat = scData_Train.CELL2CAT
        # self.n_gene = scDataset.n_vars
        # self._load_model()
        # self.model.eval()
        # emb_cells = []

        # for item in tqdm(self.scDataLoader):
        #     counts_ = item["counts"].to(self.device).squeeze()
        #     # for _idx in range(self.adata.n_obs):
        #     #     if (self.adata.X[_idx, :] == counts_.cpu().numpy()[1]).all():
        #     #         print(_idx)
        #     emb_cell = self.model.extra_repr(counts_)
        #     emb_cells.append(emb_cell.detach().cpu().numpy())

        # emb_cells = np.concatenate(emb_cells, axis=0)
        # self.adata.obsm["Islander"] = emb_cells
        
        # reducer = umap.UMAP(min_dist=0.5)
        # embedding = reducer.fit_transform(self.adata.obsm["Islander"])
        # self.adata.obsm["Islander_UMAP"] = embedding

        # self._save_adata_()
        return



if __name__ == "__main__":
    #
    args = Parser_Benchmarker()
    benchmarker = scIB(args=args)
    benchmarker._benchmark_()

    # jaxlib is version 0.4.7, but this version of jax requires version >= 0.4.14.
