from ._pystrings import _format_string_printing_font as format_pystring
import numpy as np
import os, shutil

def _enumerate_search_keyword_over_path(
    path_dir_to_search,
    keyword=None,
    search_files=False,
    bypass_list=[".anaconda3", ".local"],
):

    """
    Enumerates over all directories, sub-directories, and files within a given search directory.
    Appends results to a list and returns.

    Parameters:
    -----------
    path_dir_to_search
        type:str

    keyword
        type:str
    
    search_files
        default:False
       
    bypass_list (optional)
        type:list
        
    Returns:
    --------
    found_paths
    """

    matching_dirpath, matching_files = [], []
    for dirpath, dirnames, filenames in os.walk(path_dir_to_search):
        if keyword in dirpath:
            chunks = np.array(dirpath.split("/"))
            if np.any([chunk in bypass_list for chunk in chunks]):
                continue
            else:
                matching_dirpath.append(dirpath)
        if search_files:
            if keyword in filenames:
                matching_files.append(filenames)

    if search_files:
        return matching_dirpath, matching_files
    print(
        "{} directories were identified containing keyword: {}".format(
            format_pystring(str(len(matching_dirpath)), ["RED", "BOLD"]),
            format_pystring(keyword, ["RED", "BOLD"]),
        )
    )
    return matching_dirpath

def _find_and_destroy_pycache(path_dir_to_search="/home/mvinyard/"):

    __pycache__to_destroy = _enumerate_search_keyword_over_path(
        path_dir_to_search,
        keyword="__pycache__",
        search_files=False,
        bypass_list=[".anaconda3", ".local"],
    )
    print("Deleting...", flush=True, end="\r")
    [shutil.rmtree(pycache) for pycache in __pycache__to_destroy]
    print("Deleting... Done.")