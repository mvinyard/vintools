from ..._plotting import share_seq
from ..._plotting import vin_colors

import os
from glob import glob
import numpy as np
import anndata as a
import matplotlib.pyplot as plt
import scanpy as sc
from harmony import harmonize

import warnings

warnings.filterwarnings("ignore")


def _get_sample_name(sample_path):

    sample_name = os.path.basename(sample_path).split("_")[0]

    return sample_name


def _scanpy_filter(adata, n_genes_cutoff=2500, pct_mt_cutoff=5):

    adata = adata[adata.obs.n_genes_by_counts < n_genes_cutoff, :]
    adata = adata[adata.obs.pct_counts_mt < pct_mt_cutoff, :]


def _run_scanpy_standard_qc(adata, plot=True):

    """"""
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    if plot:

        sc.pl.highest_expr_genes(
            adata, n_top=20,
        )
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
        )

        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt")
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")


def _scanpy_normalize(
    adata, plot=True, perform_filtering=False, regress_out=False, max_value=10
):

    print("Running scanpy normalization against counts total...")
    sc.pp.normalize_total(adata, target_sum=1e4)
    print("Running scanpy log normalization...")
    sc.pp.log1p(adata)

    print("Identifying highly variable genes...")
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    if plot:
        sc.pl.highly_variable_genes(adata)

    if perform_filtering or regress_out:
        adata.raw = adata

        if perform_filtering:
            adata = adata[:, adata.var.highly_variable]

        if regress_out:
            sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

    print("Scaling gene expression...")
    sc.pp.scale(adata, max_value=max_value)


def _run_scanpy_pca_umap(
    adata, svd_solver="arpack", plot_report=False, plot_genes=False
):

    """

    plot_genes: list
    """

    # PCA
    print("Running Scanpy PCA...")
    sc.tl.pca(adata, svd_solver=svd_solver)
    if plot_genes:
        sc.pl.pca(adata, color=plot_genes)
    if plot_report:
        sc.pl.pca_variance_ratio(adata, log=True)

    # UMAP
    print("Running Scanpy UMAP...")
    sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
    sc.tl.umap(adata)
    if plot_genes:
        # just a quick check
        sc.pl.umap(adata, color=plot_genes)


def _run_clustering(adata, leiden_resolution=1, louvain_resolution=1, plot=False):

    print("Performing Leiden clustering (resolution={})...".format(leiden_resolution))
    sc.tl.leiden(adata, resolution=leiden_resolution)
    print("Performing Louvain clustering (resolution={})...".format(louvain_resolution))
    sc.tl.louvain(adata, resolution=louvain_resolution)
    if plot:
        assert type(plot) == list, print("plot must be list")
        sc.pl.umap(adata, color=plot)
    else:
        sc.pl.umap(adata, color=["leiden", "louvain", "sample"])


def _run_Harmony(adata, harmonize_on="X_pca", batch_key="sample"):

    """"""

    if harmonize_on == "X_pca":
        print(
            "running harmony on adata.obsm['X_pca']; shape: {}".format(
                adata.obsm["X_pca"].shape
            )
        )
        adata.obsm["X_harmony_on_X_pca"] = harmonize(
            adata.obsm["X_pca"], adata.obs, batch_key=batch_key
        )
    if harmonize_on == "X":
        print(
            "running harmony on adata.X; shape: {}\n\n{}: Time and memory requirements will be intensified when harmonizing across adata.X compared to a dimensional reduction.".format(
                adata.X.shape, v.ut.format_pystring("WARNING", ["RED", "BOLD"])
            )
        )
        adata.obsm["X_harmony_on_X"] = harmonize(
            adata.X, adata.obs, batch_key=batch_key
        )


def _run_umap(adata, obsm_key="X_harmony_on_X_pca"):

    """"""
    import umap
    from sklearn.preprocessing import StandardScaler

    print("Calculating umap embedding...")

    reducer = umap.UMAP()
    scaled_data = StandardScaler().fit_transform(adata.obsm[obsm_key])
    adata.obsm["X_umap_harmony"] = embedding = reducer.fit_transform(scaled_data)


def _singlefig_umap_presets(
    title="UMAP",
    figsize=(8, 6.5),
    title_fontsize=12,
    label_fontsize=10,
    title_adjust_y=1.05,
):

    """"""

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1)

    ax.set_title(title, fontsize=16, y=1.05)
    ax.set_xlabel("UMAP-1", fontsize=label_fontsize)
    ax.set_ylabel("UMAP-2", fontsize=label_fontsize)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    ax.axes.xaxis.set_ticks([])
    ax.axes.yaxis.set_ticks([])

    return fig, ax


def _quick_extend_color_palette(n_ids):

    extended_default = share_seq()["colors"].tolist()
    if n_ids > len(extended_default):
        for color in range(n_ids - len(extended_default)):
            extended_default.append(np.tile(vin_colors(separable=True), 3)[color])

        return extended_default


def _plot_UMAP_by_categorical_obskey(
    adata, obs_key=None, legend_cols=2, separate_colors=False
):

    """"""

    if not adata.obs.index.values.dtype == "int64":
        adata.obs = adata.obs.reset_index()

    assert obs_key != None, print("AssertionError: Define obskey !")

    n_ids = len(adata.obs[obs_key].unique())

    if int(n_ids) > share_seq()["colors"].shape[0]:
        colors = _quick_extend_color_palette(n_ids)
    else:
        colors = share_seq()["colors"]

    embedding = adata.obsm["X_umap_harmony"]

    try:
        adata.obs[obs_key] = adata.obs[obs_key].astype(int)
    except:
        pass

    for n, sample in enumerate(np.sort(adata.obs[obs_key].unique())):
        sample_idx = adata.obs.loc[adata.obs[obs_key] == sample].index.astype(int)
        emb_x, emb_y = embedding[:, 0][sample_idx], embedding[:, 1][sample_idx]
        if separate_colors:
            n *= 4
        plt.scatter(emb_x, emb_y, label=sample, s=2, c=colors[n], alpha=1)

    plt.legend(
        markerscale=3,
        edgecolor="w",
        fontsize=14,
        ncol=legend_cols,
        handletextpad=None,
        bbox_to_anchor=(0.5, 0.0, 0.95, 1),
    )


def _plot_UMAP_by_GEX(adata, gene=None, cmap="Purples"):

    """"""

    assert gene != None, print("AssertionError: Define gene !")

    gene_vals = adata[:, gene].X
    embedding = adata.obsm["X_umap_harmony"]

    plt.scatter(
        embedding[:, 0], embedding[:, 1], c="black", cmap=cmap, s=8, alpha=0.5, zorder=0
    )
    plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c="whitesmoke",
        cmap=cmap,
        s=4,
        alpha=0.5,
        zorder=0,
    )

    plt.scatter(
        embedding[:, 0],
        embedding[:, 1],
        c=gene_vals,
        cmap=cmap,
        s=5,
        alpha=0.5,
        zorder=2,
    )
    plt.colorbar(shrink=0.35, aspect=10)


def _adjust_clustering_resolution(
    adata,
    algorithm="both",
    desired_resolution=None,
    desired_adjustment=1,
    title_fontsize=12,
    label_fontsize=10,
    figsize=(8, 6.5),
    title_adjust_y=1.05,
    legend_cols=2,
    separate_colors=False,
):

    """
    adata

    algorithm
        louvain, leiden, both
        default: both

    desired_resolution
        Min: 0, Max: 1
        default: None

    desired_adjustment
        default: 1 (no adjustment)
    """

    current_resolution = float(adata.uns["louvain"]["params"]["resolution"])

    if desired_resolution is None:
        print("Adjusting resolution by a factor of {}".format(desired_adjustment))
        desired_resolution = current_resolution * desired_adjustment
        print("Adjusted resolution: {}".format(desired_resolution))

    if algorithm == "leiden" or algorithm == "both":
        print(
            "Performing Leiden clustering at resolution {}...".format(
                desired_resolution
            )
        )
        sc.tl.leiden(adata, resolution=desired_resolution)

        plot_title = "Harmony Integration: Leiden\nadjusted resolution: {}".format(
            desired_resolution
        )
        _singlefig_umap_presets(
            title=plot_title,
            figsize=figsize,
            title_fontsize=title_fontsize,
            label_fontsize=label_fontsize,
            title_adjust_y=title_adjust_y,
        )

        _plot_UMAP_by_categorical_obskey(
            adata,
            obs_key="leiden",
            legend_cols=legend_cols,
            separate_colors=separate_colors,
        )

    if algorithm == "louvain" or algorithm == "both":
        print(
            "Performing Louvain clustering at resolution {}...".format(
                desired_resolution
            )
        )
        sc.tl.louvain(adata, resolution=desired_resolution)
        plot_title = "Harmony Integration: Louvain\nadjusted resolution: {}".format(
            desired_resolution
        )
        _singlefig_umap_presets(
            title=plot_title,
            figsize=figsize,
            title_fontsize=title_fontsize,
            label_fontsize=label_fontsize,
            title_adjust_y=title_adjust_y,
        )
        _plot_UMAP_by_categorical_obskey(
            adata,
            obs_key="louvain",
            legend_cols=legend_cols,
            separate_colors=separate_colors,
        )
    print(
        "Previous resolution: {}\nAdjusted resolution: {}".format(
            current_resolution, desired_resolution
        )
    )


def _scanpy_dotplots_GEX(adata, marker_genes, clustering_algorithm="leiden"):

    """"""

    adata.obs[clustering_algorithm] = adata.obs[clustering_algorithm].astype(str)
    adata.obs[clustering_algorithm] = adata.obs[clustering_algorithm].astype(str)

    sc.pl.dotplot(adata, marker_genes, groupby=clustering_algorithm)
