# utilities

from ._vcf_tools import read_vcf
from ._send_email import _send_email as email
from ._cell_complete_audio import _ipy_cell_complete_audio as done
from ._flatten_lists import _flatten_list_of_lists as flatten
# from ._zero_centered_random import _zero_centered_random as zero_centered_random
from ._pystrings import _format_string_printing_font as format_pystring
from ._search_and_destroy import _find_and_destroy_pycache as findNdestroy_pycache
from ._pyGSUTILS_module._pyGSUTILS import _pyGSUTILS as pyGSUTILS

from ._dynamical_import_of_function_from_string import _dynamical_import_of_function_from_string as import_from_string
from ._deinterleaf import _deinterleaf_fastq as deinterleaf_fastq
from ._report_imports import _report_imports as imports
from ._get_vintools_version import _get_vintools_version as version

from ._smooth_data import _partition as partition
from ._smooth_data import _smooth as smooth

from ._flexible_mkdir import _flexible_mkdir as mkdir_flex


from ._get_data_bounds import _get_data_bounds as get_data_bounds

# from ._time import _time as time
from ._fetch_cpu_count import _fetch_cpu_count as fetch_n_cpus
from ._use_n_cores import _use_n_cores as use_n_cores


from ._torch_device import _set_device as set_device
from ._torch_device import _torch_device as torch_device