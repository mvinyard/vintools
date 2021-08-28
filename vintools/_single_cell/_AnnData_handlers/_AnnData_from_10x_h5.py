import h5py
import anndata
import numpy as np
import pandas as pd
import scipy.sparse

def _get_var_df(mat_group, var_col=["id", "genome", "feature_type"]):

    var_df = pd.DataFrame()

    for col in var_col:
        var_df[col] = np.array(mat_group["features"][col]).astype(str)
    return var_df


def _get_matrix_from_h5(mat_group):

    data = np.array(mat_group["data"])
    indices = np.array(mat_group["indices"])
    indptr = np.array(mat_group["indptr"])
    shape = np.array(mat_group["shape"])
    matrix = scipy.sparse.csc_matrix((data, indices, indptr), shape=shape).T

    return matrix


def _get_barcodes(mat_group):

    return pd.DataFrame(
        mat_group["barcodes"][()].astype(str), columns=["cell_barcodes"]
    )

def _AnnData_from_10x_h5(
    h5_path, var_col=["id", "genome", "feature_type"], silent=False
):

    """

    Parameters:
    -----------
    h5_file

    Returns:
    --------
    adata
        type: anndata._core.anndata.AnnData

    Notes:
    ------
    (1)  Assumes h5_file['matrix'] exists. Current output of 10x CellRanger.

    (2)  Developed with input from the following sources:
            https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
            https://stackoverflow.com/questions/40388792/how-to-decode-a-numpy-array-of-encoded-literals-strings-in-python3-attributeerr
            https://stackoverflow.com/questions/28170623/how-to-read-hdf5-files-in-python/41586571
    """

    h5_file = h5py.File(h5_path, mode="r")

    adata = anndata.AnnData(_get_matrix_from_h5(h5_file["matrix"]))
    adata.var = _get_var_df(h5_file["matrix"], var_col=var_col)
    adata.obs = _get_barcodes(h5_file["matrix"])

    if not silent:
        print(adata)

    return adata