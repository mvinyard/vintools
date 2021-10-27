# utilities

from ._data_utils._create_empty_dict import _create_EmptyDict as EmptyDict

from ._data_utils._write_multi_df_to_excel import _write_multi_df_to_excel as df_to_excel
from ._data_utils._read_txt import _read_txt as read_txt
from ._data_utils._vcf_tools import read_vcf
from ._data_utils._flatten_lists import _flatten_list_of_lists as flatten
from ._data_utils._get_data_bounds import _get_data_bounds as get_data_bounds
from ._data_utils._torch_device import _set_device as set_device
from ._data_utils._torch_device import _torch_device as torch_device
from ._data_utils._pyGSUTILS_module._pyGSUTILS import _pyGSUTILS as pyGSUTILS
from ._data_utils._deinterleaf import _deinterleaf_fastq as deinterleaf_fastq
from ._data_utils._smooth_data import _partition as partition
from ._data_utils._smooth_data import _smooth as smooth
from ._data_utils._parse_gtf import _parse_gtf_as_dictionary as parse_gtf
from ._data_utils._FileHandler import _FileHandler as FileHandler
from ._data_utils._check_fix_file_extension import _check_fix_file_extension as secure_file_extension
from ._data_utils._glob_dict import _glob_dict as glob_dict
from ._data_utils._fix_path_for_glob import _fix_path_for_glob as secure_glob_path

# from ._zero_centered_random import _zero_centered_random as zero_centered_random
from ._ux_utils._pystrings import _format_string_printing_font as format_pystring
from ._ux_utils._cell_complete_audio import _ipy_cell_complete_audio as done
from ._ux_utils._print_underline import _print_underline as print_underlined

from ._system_utils._get_basename_no_extension import  _get_basename_no_extension as basename_no_ext
from ._system_utils._report_imports import _report_imports as imports
from ._system_utils._get_vintools_version import _get_vintools_version as version
from ._system_utils._fetch_cpu_count import _fetch_cpu_count as fetch_n_cpus
from ._system_utils._use_n_cores import _use_n_cores as use_n_cores
from ._system_utils._clone_GitHub_repo import _clone_GitHub_repo as git_clone
from ._system_utils._get_pypi_package_loc import _get_pypi_package_loc as which
from ._system_utils._flexible_mkdir import _flexible_multilevel_mkdir as mkdir_flex
from ._system_utils._search_and_destroy import (
    _find_and_destroy_pycache as findNdestroy_pycache,
)
from ._system_utils._send_email import _send_email as email
from ._system_utils._dynamical_import_of_function_from_string import (
    _dynamical_import_of_function_from_string as import_from_string,
)


# from ._time import _time as time
