from ._fetch_cpu_count import _fetch_cpu_count


def _use_n_cores(n_cores):

    """
    Given a number of desired cores to be used for a process, this function checks how many cores are available and 
    returns the whatever value is greater, cores available or cores requested. 
    
    Parameters:
    -----------
    n_cores
        number of cores desired to be used for a given process
        
    Returns:
    --------
    n_cores
        number of cores to be used for a given process
        type: int
    
    Notes:
    ------
    """
    n_cpus_available = _fetch_cpu_count()

    if n_cores > n_cpus_available:
        print(
            "{} cores requested. Only {} available! Using {}".format(
                n_cores, n_cpus_available, n_cpus_available
            )
        )
        n_cores = n_cpus_available
    else:
        print("Using {} cores...".format(n_cores))

    return n_cores
