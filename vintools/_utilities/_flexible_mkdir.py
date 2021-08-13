
# package imports #
# --------------- #
import os


def _flexible_mkdir(path):
    
    """
    Create a directory or ignore if already present. 
    
    Parameters:
    -----------
    path
    
    Returns:
    --------
    None
        os.mkdir(path) (or nothing)
    
    Notes:
    ------
    """
    
    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)