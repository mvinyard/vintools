# utilities

from ._data_utils._vcf_tools import read_vcf
from ._system_utils._send_email import _send_email as email
from ._ux_utils._cell_complete_audio import _ipy_cell_complete_audio as done
from ._data_utils._flatten_lists import _flatten_list_of_lists as flatten
# from ._zero_centered_random import _zero_centered_random as zero_centered_random
from ._ux_utils._pystrings import _format_string_printing_font as format_pystring
from ._system_utils._search_and_destroy import _find_and_destroy_pycache as findNdestroy_pycache
from ._data_utils._pyGSUTILS_module._pyGSUTILS import _pyGSUTILS as pyGSUTILS

from ._system_utils._dynamical_import_of_function_from_string import _dynamical_import_of_function_from_string as import_from_string
from ._data_utils._deinterleaf import _deinterleaf_fastq as deinterleaf_fastq
from ._system_utils._report_imports import _report_imports as imports
from ._system_utils._get_vintools_version import _get_vintools_version as version

from ._data_utils._smooth_data import _partition as partition
from ._data_utils._smooth_data import _smooth as smooth

from ._system_utils._flexible_mkdir import _flexible_mkdir as mkdir_flex


from ._data_utils._get_data_bounds import _get_data_bounds as get_data_bounds

# from ._time import _time as time
from ._system_utils._fetch_cpu_count import _fetch_cpu_count as fetch_n_cpus
from ._system_utils._use_n_cores import _use_n_cores as use_n_cores


from ._data_utils._torch_device import _set_device as set_device
from ._data_utils._torch_device import _torch_device as torch_device


from ._system_utils._clone_GitHub_repo import _clone_GitHub_repo as git_clone
from ._system_utils._get_pypi_package_loc import _get_pypi_package_loc as which

from ._data_utils._parse_gtf import _parse_gtf_as_dictionary as parse_gtf