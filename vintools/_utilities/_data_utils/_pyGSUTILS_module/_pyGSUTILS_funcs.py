# gsutil_funcs.py
import subprocess
import os, glob
from ..._ux_utils._pystrings import _format_string_printing_font

def _init_funcs():

    """"""

    base_gsutil_command = "gsutil"
    cwd = os.getcwd()
    gsutil_loc = (
        subprocess.check_output("which gsutil", shell=True).decode().strip("\n")
    )

    return base_gsutil_command, cwd, gsutil_loc

def _get_local_copied_file_path(infile_path, download_dir):
    
    """"""
    
    filename = os.path.basename(infile_path)
    local_filepath = os.path.join(download_dir, infile_path)
    
    return local_filepath

def _decode_bytes_to_string_list(bytes_string):

    decoded_string_list = []

    for string in bytes_string.split(b"\n"):
        decoded_string_list.append(string.decode())

    return decoded_string_list


def _list_available_buckets(silent=False, return_buckets=True):

    """
    Returns a list of buckets from which data can be retrieved.
    Default behavior prints list. Can be returned as an object.

    Parameters:
    -----------
    silent
        default: False

    return_buckets
        default: False

    Returns:
    --------
    [optional] buckets
    """

    base_ls_command = "gsutil ls gs://"
    bucket_bytes = subprocess.check_output(base_ls_command, shell=True)
    buckets = _decode_bytes_to_string_list(bucket_bytes)

    if not silent:
        for bucket in buckets:
            print(bucket)

    if return_buckets:
        return buckets


def _read_cloud_storage_dir(path=None, silent=False, return_ls=True):

    """
    By default, lists all project buckets.
    """

    base_ls_command = "gsutil ls gs://"

    if path != None:
        command_path = os.path.join(base_ls_command, (path + "*"))
    else:
        command_path = base_ls_command

    ls_read_data_bytes = subprocess.check_output(command_path, shell=True)
    ls_outs = _decode_bytes_to_string_list(ls_read_data_bytes)

    if not silent:
        for item in ls_outs:
            print(item)

    if return_ls:
        return ls_outs
    
def _gsutil_cp(base_command, source_path=None, destination_path='./', multithread=True, recursive=True, verbose=False):
    
    """"""
    
    assert source_path != None, print("Define the path to an object for download!")
    
    if multithread:
        cp_command = " ".join([base_command, "-m cp"])
    else:
        cp_command = " ".join([base_command, "cp"])
        
    if recursive:
        cp_command = " ".join([cp_command, "-r"])
        
    cp_command = "".join([cp_command, " gs://", source_path])
    destination_path = os.path.join(destination_path, os.path.basename(source_path))
    cp_command = " ".join([cp_command, destination_path])
    
    if verbose:
        print(_format_string_printing_font("gsutil cp command:", ["BOLD", "CYAN"]), cp_command)
    os.system(cp_command)