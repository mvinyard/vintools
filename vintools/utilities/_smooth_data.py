
# package imports #
# --------------- #
import numpy as np

def _partition(unpartitioned_items, groupsize=5):
    
    """
    Yield successive n-sized chunks. 

    Parameters:
    -----------
    unpartitioned_items

    groupsize

    Returns:
    --------
    list of np.arrays of size N (except perhaps the last array)

    Notes:
    ------
    source: https://www.geeksforgeeks.org/break-list-chunks-size-n-python/
    """

    for i in range(0, len(unpartitioned_items), groupsize):
        yield np.array(unpartitioned_items[i : i + groupsize])


def _smooth(unpartitioned_items, groupsize=5):
    
    """
    Given a list and groupsize, returns smoothed mean and stdev over groupsize-partitioned list. 

    Parameters:
    -----------
    unpartitioned_items

    groupsize

    Returns:
    --------
    smoothed_means
        type: np.ndarray

    smoothed_stdev
        type: np.ndarray

    Notes:
    ------
    """
        
    smoothed_means = []
    smoothed_stdev = []

    partitioned_list = list(_partition(unpartitioned_items, groupsize))
    for array in partitioned_list:
        smoothed_means.append(array.mean())
        smoothed_stdev.append(array.std())

    return np.array(smoothed_means), np.array(smoothed_stdev)