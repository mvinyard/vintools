import os
import shutil


def _create_data_cache(data_cache_dir="./.data_cache/", project=None):
    
    """
    Given a path, makes a local cache directory for downloaded data. 
    """
    
    if not os.path.exists(data_cache_dir):
        if project != None:
            project = "for {}...".format(project)
        
        print("Creating data cache directory at {} {}".format(data_cache_dir, project))
        os.mkdir(data_cache_dir)
            
    return data_cache_dir

def _remove_cache_dirs(path_dir_to_search="./", cache_id_keyword="data_cache"):
    
    cache_dir_to_destroy = _enumerate_search_keyword_over_path(
        path_dir_to_search,
        keyword=cache_id_keyword,
        search_files=False,
        bypass_list=[".anaconda3", ".local"],
    )
    print("Deleting...", flush=True, end="\r")
    [shutil.rmtree(data_cache) for data_cache in cache_dir_to_destroy]
    print("Deleting... Done.")