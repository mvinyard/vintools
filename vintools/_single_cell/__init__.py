# single-cell data utilities __init__.py

__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# general data formatting functions
from ._general_formatting_tools import (
    _check_if_scipy_sparse_mtx as check_if_scipy_sparse_mtx,
)

# single-cell AnnData formatting functions
from ._AnnData_handlers._format_AnnData import _format_AnnData as format_adata
from ._Scorch_Module._scorch import Scorch
from ._AnnData_handlers._AnnData_from_10x_h5 import _AnnData_from_10x_h5 as read_10x_h5

# 10x cellranger auxiliary functions
from ._10x_cellranger_auxiliary._post_10x_pipestance import (
    _post_10x_pipestance as post_10x_pipestance,
)
from ._10x_cellranger_auxiliary._post_10x_pipestance import (
    _run_10x_post_pipestance as run_10x_post_pipestance,
)
