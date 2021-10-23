
# package imports #
# --------------- #
import os

# local imports #
# ------------- #
from .._ux_utils._pystrings import _format_string_printing_font


def _flexible_mkdir(path, verbose):

    if os.path.exists(path):
        pass
    else:
        os.mkdir(path)
        if verbose:
            msg = _format_string_printing_font("Directory created", ["BOLD", "CYAN"])
            print("{}: {}".format(msg, path))


def _flexible_multilevel_mkdir(path, verbose=False):

    """
    Create a directory or ignore if already present. Can create multiple levels of directories.

    Parameters:
    -----------
    path

    Returns:
    --------
    None
        os.mkdir(path) (or nothing)

    Notes:
    ------
    (1) supports multi-level directory creation

    """

    parsed_path = []

    for n, directory in enumerate(path.split("/")):        
        parsed_path.append(directory)
        if len(parsed_path) == len(path.split("/")):
            _flexible_mkdir("/".join(parsed_path), verbose)
        else:
            _flexible_mkdir("/".join(parsed_path), verbose=False)
            