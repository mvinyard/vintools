import multiprocessing

def _fetch_cpu_count():
    
    """
    Returns the number of available CPUs on machine in use. 
    
    Parameters:
    -----------
    None
    
    Returns:
    --------
    multiprocessing.cpu_count()
    
    Notes:
    ------
    None
    """

    return multiprocessing.cpu_count()
