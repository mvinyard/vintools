
import numpy as np
import os

from ..._utilities._flexible_mkdir import _flexible_mkdir
from ..._utilities._get_pypi_package_loc import _get_pypi_package_loc

def _choose_n_cancer_samples(df_meta, n_samples):

    cancer_samples = df_meta.loc[df_meta.condition == "cancer"]["sample"].values
    return np.random.choice(cancer_samples, n_samples, replace=False)


def _subset_input_dfs(df, df_meta, n_cancer_samples):

    """
    Parameters:
    -----------

    Returns:
    --------

    Notes:
    ------
    """

    normal_samples = df_meta.loc[df_meta.condition == "normal"]["sample"].values
    cancer_samples = _choose_n_cancer_samples(df_meta, n_cancer_samples)

    df_meta_subset = df_meta.loc[
        df_meta["sample"].isin(np.append(normal_samples, cancer_samples))
    ]

    df_subset = df[df_meta_subset["sample"].values]

    return df_subset, df_meta_subset


def _prep_DESeq2_input_data(dfe, df_meta, n_cancer_samples, tmp_run_dir):

    """"""

    df_, dfm_ = _subset_input_dfs(dfe, df_meta, n_cancer_samples)
    
    path_data = os.path.join(tmp_run_dir, "gex_counts.csv")
    path_meta = os.path.join(tmp_run_dir, "gex_meta.csv")

    df_.to_csv(path_data)
    dfm_.to_csv(path_meta, index_label=False)


def _run_DESeq2(
    dfe, df_meta, n_cancer_samples, outs_dir="./", run_name="testRun"
):

    """"""
    
    software_dir = os.path.dirname(_get_pypi_package_loc())
    path = os.path.join(software_dir, "vintools/vintools/_tools/_DESeq2/DESeq2.R")
    
    tmp_run_dir="tmp_{}".format("".join(np.random.choice(list(string.ascii_lowercase), 6)))
    
    _flexible_mkdir(outs_dir)
    _flexible_mkdir(tmp_run_dir)
    _prep_DESeq2_input_data(dfe, df_meta, n_cancer_samples, tmp_run_dir)

    executable = " ".join(["Rscript", path, outs_dir, run_name])
    print(executable)
    os.system(executable)
    
    os.system("rm -r {}".format(tmp_run_dir))


def _run_batch_DESeq2(dfe, df_meta, n_iters, n_cancer_samples):
    for iteration in range(n_iters):
        outsdirpath = "nsamples_{}".format(str(n_cancer_samples))
        _run_DESeq2(
            dfe,
            df_meta,
            n_cancer_samples,
            outs_dir=outsdirpath,
            run_name="iter_{}".format(str(iteration)),
        )
