
# single-cell data utilities __init__.py

__author__ = ', '.join([
    'Michael E. Vinyard'
])
__email__ = ', '.join([
    'vinyard@g.harvard.edu',
])

# general data formatting functions
from ._general_formatting_tools import _check_if_scipy_sparse_mtx as check_if_scipy_sparse_mtx

# single-cell AnnData formatting functions
from ._AnnData_handlers._format_AnnData import _format_AnnData as format_adata
from ._Scorch_Module._scorch import Scorch
