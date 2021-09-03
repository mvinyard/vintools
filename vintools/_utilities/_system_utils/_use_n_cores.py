
from ._fetch_cpu_count import _fetch_cpu_count

def _use_n_cores(n_cores):
    
    """"""
    n_cpus_available = _fetch_cpu_count()
    
    if n_cores > n_cpus_available:
        print("{} cores requested. Only {} available! Using {}".format(n_cores, n_cpus_available, n_cpus_available))
        n_cores = n_cpus_available
    else:
        print("Using {} cores...".format(n_cores))
