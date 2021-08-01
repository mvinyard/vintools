
from ._scorch_supporting_functions import _get_sample_name, _scanpy_filter, _run_scanpy_standard_qc, _scanpy_normalize, _run_scanpy_pca_umap, _run_clustering, _run_Harmony, _run_umap, _singlefig_umap_presets, _quick_extend_color_palette, _plot_UMAP_by_categorical_obskey, _plot_UMAP_by_GEX, _adjust_clustering_resolution, _scanpy_dotplots_GEX

import time
import glob
import scanpy as sc
import anndata as a
    
class Scorch:

    def __init__(self):

        """"""        
        
        import harmony
        self.start_time = time.time()

        print(
            "Harmony version: {}\nScanpy version: {}\n".format(
                harmony.__version__, sc.__version__
            )
        )

    def read_plural_10x_h5(self, data_dir, return_adata=False):
    
        sample_names, adata_list = [], []
        h5_paths = glob.glob(data_dir)

        count = 1
        for path in h5_paths:
            sample_name = _get_sample_name(path)
            sample_names.append(sample_name)
            sample_adata = sc.read_10x_h5(path)
            sample_adata.var_names_make_unique()
            sample_adata.obs.index = sample_adata.obs.index.str.replace(
                "-1", "-{}".format(sample_name)
            )
            adata_list.append(sample_adata)
            count += 1
        self.adata = a.concat(
            adata_list, merge="same", label="sample", keys=sample_names
        )
        if return_adata:
            return self.adata

    def scanpy_qc(
        self, plot=True, perform_filter=False, n_genes_cutoff=2500, pct_mt_cutoff=5
    ):

        """"""

        _run_scanpy_standard_qc(self.adata, plot=plot)

        if perform_filter:
            print(
                "Filtering AnnData with a minimum gene count cutoff of: {} and minimum %MT reads of: {}%".format(
                    n_genes_cutoff, pct_mt_cutoff
                )
            )
            _scanpy_filter(
                self.adata, n_genes_cutoff=n_genes_cutoff, pct_mt_cutoff=pct_mt_cutoff
            )

    def scanpy_counts_analysis(
        self, plot=True, perform_filtering=False, regress_out=False
    ):

        """"""

        _scanpy_normalize(
            self.adata,
            plot=plot,
            perform_filtering=perform_filtering,
            regress_out=regress_out,
        )

    def scanpy_dimensional_reduction_clustering(
        self,
        svd_solver="arpack",
        plot_pca_report=False,
        plot_genes=False,
        leiden_resolution=1,
        louvain_resolution=1,
        plot=False,
    ):

        _run_scanpy_pca_umap(
            self.adata,
            svd_solver=svd_solver,
            plot_report=plot_pca_report,
            plot_genes=plot_genes,
        )
        _run_clustering(
            self.adata,
            leiden_resolution=leiden_resolution,
            louvain_resolution=louvain_resolution,
            plot=plot,
        )

    def harmonize(self, harmonize_on="X_pca", batch_key="sample"):
        _run_Harmony(self.adata, harmonize_on="X_pca", batch_key="sample")

    def umap(self, obsm_key="X_harmony_on_X_pca"):
        _run_umap(self.adata, obsm_key=obsm_key)

    def visualize(
        self,
        plot_by=None,
        title=None,
        title_fontsize=12,
        label_fontsize=10,
        cmap="Purples",
        figsize=(8, 6.5),
        title_adjust_y=1.05,
        legend_cols=2,
        separate_colors=False,
    ):

        """"""

        if plot_by is None:
            print("plot-type:None")
        elif plot_by in self.adata.obs.columns:
            plot_type = "categorical"
        elif plot_by in self.adata.var_names:
            plot_type = "gene"
        else:
            print("Plot-type not well defined or perhaps gene name is inaccurate.")

        if title is None:
            title = plot_by
        _singlefig_umap_presets(title=title, figsize=figsize, title_fontsize=title_fontsize, label_fontsize=label_fontsize, title_adjust_y=title_adjust_y,)

        if plot_type == "categorical":
            _plot_UMAP_by_categorical_obskey(self.adata, obs_key=plot_by, legend_cols=legend_cols, separate_colors=separate_colors)

        if plot_type == "gene":
            _plot_UMAP_by_GEX(self.adata, gene=plot_by, cmap=cmap)

    def adjust_clustering_resolution(self, algorithm="both", 
                                     desired_resolution=None, 
                                     desired_adjustment=1, 
                                     title_fontsize=12, 
                                     label_fontsize=10, 
                                     figsize=(8, 6.5), 
                                     title_adjust_y=1.05, 
                                     legend_cols=2, 
                                     separate_colors=False):

        _adjust_clustering_resolution(self.adata, algorithm=algorithm, 
                                      desired_resolution=desired_resolution, 
                                      desired_adjustment=1, 
                                      title_fontsize=title_fontsize, 
                                      label_fontsize=label_fontsize, 
                                      figsize=figsize, 
                                      title_adjust_y=title_adjust_y, 
                                      legend_cols=legend_cols, 
                                      separate_colors=separate_colors)

    def cluster_dotplot(self, marker_genes=None, clustering_algorith="leiden"):
        
        if marker_genes is None:
            marker_genes = self.marker_genes
            print("Using stored marker genes...")
        _scanpy_dotplots_GEX(self.adata, marker_genes, clustering_algorithm=clustering_algorith)
        
        
    def check_time(self):
    
        self.end_time = time.time()
        time_elapsed = (self.end_time - self.start_time) / 60
        print("{:.2f} mins.".format(time_elapsed))
        
