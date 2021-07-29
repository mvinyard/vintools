# from ._general_formatting_tools import _check_if_scipy_sparse_mtx

def _format_AnnData_mtx_as_numpy_array(adata):

    """
    Checks:
        (1) scipy.stats.sparse_matrix format
        (2) numpy.ndarray format
    Then:
        Converts adata.X to numpy.ndarray or reports that the array was not able to be converted. 
    
    Parameters:
    -----------
    adata

    Returns:
    --------
    None
        AnnData is modified in place. adata.X is converted to a numpy.ndarray.
    """
    
    import vintools as v
    
    array_is_not_status_list = ["\nadata.X is not:"]
    
    # check if adata.X is a scipy sparse matrix of some sort    
    scipy_status = v.sc.check_if_scipy_sparse_mtx(adata.X)
    if scipy_status in ["coo","csr"]:
        out_message="scipy.sparse.{}_matrix".format(scipy_status)
        print("adata.X was of dtype: {}.".format(v.ut.format_pystring(out_message,["RED","BOLD"] )))
        adata.X = adata.X.toarray()
        print("adata.X has been converted to dtype: {}.".format(v.ut.format_pystring("numpy.ndarray", ["RED","BOLD"])))
        return None
    else:
        array_is_not_status_list.append("scipy.sparse.csr_matrix")
        array_is_not_status_list.append("scipy.sparse.coo_matrix")
        
    # check if numpy array
    if type(adata.X).__name__ == 'ndarray':
        out_message="adata.X"
        print("{} is already of dtype: numpy.ndarray".format(v.ut.format_pystring(out_message, ["RED","BOLD"])))
    else:
        try:
            adata.X = adata.X.toarray()
            out_message="adata.X"
            print("{} converted to dtype: numpy.ndarray".format(v.ut.format_pystring(out_message, ["RED","BOLD"])))
        except:
            array_is_not_status_list.append("np.ndarray")
            print("adata.X was unable to be converted to an array.")
            for i, term in enumerate(array_is_not_status_list):
                if i != 0:
                    print("\n\t{}".format(v.ut.format_pystring(term, ["BOLD", "RED"])))
                else:
                    print(term)
                    
                    
    class _format_AnnData:

        """
        Parameters:
        -----------
        adata

        Returns:
        --------
        None
            AnnData object is modified in place.
        """

        def __init__(self, adata):

            self.adata = adata

        def X_to_array(self):

            _check_format_AnnData_mtx(self.adata)