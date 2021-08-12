
from ._pyGSUTILS_funcs import _init_funcs, _list_available_buckets, _read_cloud_storage_dir, _gsutil_cp, _get_local_copied_file_path

class _pyGSUTILS:
    def __init__(self):

        """
        """
        self.base_gsutil_command, self.cwd, self.gsutil_loc = _init_funcs()
        self.cp_file_paths = []
        
    def list_available_buckets(self, silent=False, return_buckets=True):

        self.available_buckets = _list_available_buckets(
            silent=silent, return_buckets=return_buckets
        )

    def read_cloud(self, path=None, silent=False, return_ls=True):

        self.ls_outs = _read_cloud_storage_dir(path=None, silent=False, return_ls=True)

    def cp(
        self, source_path=None, destination_path="./", multithread=True, recursive=True, verbose=False
    ):

        """"""

        _gsutil_cp(
            base_command=self.base_gsutil_command,
            source_path=source_path,
            destination_path=destination_path,
            multithread=multithread,
            recursive=recursive,
            verbose=verbose,
        )
        self.cp_file_paths.append(_get_local_copied_file_path(infile_path=source_path, download_dir=destination_path))