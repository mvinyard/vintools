
"""This module functions are used for general formatting of arrays."""

def _check_if_scipy_sparse_mtx(array):

    """
    Checks if input array is 'coo' or 'csr' sparse matrix. Outputs string indicating status.
    
    Parameters:
    -----------
    array
        dtype is ambiguous before running this function. 

    Returns:
    --------
    scipy.sparse type information
        ['csr','coo','not_scipy.sparse']
    """

    try:
        if array.getformat() == "csr":
            return "csr"

        if array.getformat() == "coo":
            return "coo"

    except:
        return "not_scipy.sparse"