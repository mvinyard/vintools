
# _glob_dict.py
__module_name__ = "_glob_dict.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# package imports #
# --------------- #
import glob

# local imports #
# ------------- #
from .._system_utils._get_basename_no_extension import _get_basename_no_extension
from ._fix_path_for_glob import _fix_path_for_glob
from ._FileHandler import _FileHandler


def _glob_dict(path, verbose=False):

    """
    Loads filepaths from a glob'd directory.

    Parameters:
    -----------
    path
        path to directory containing files.
        type: str

    verbose
        If true, prints FileHandler messages.
        default: False
        type: bool

    Returns:
    --------
    Dict
        keys: filenames
        values: files loaded into memory

    Notes:
    ------
    (1) Can pass with or without "/*"
    """

    filepaths = glob.glob(_fix_path_for_glob(path))

    Dict = {}
    Dict["paths"] = []

    for path in filepaths:
        Dict["paths"].append(path)
        name = _get_basename_no_extension(path)
        file = _FileHandler(filepath=path, verbose=verbose)
        Dict[name] = file.read(return_file=True)

    return Dict