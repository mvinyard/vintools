# main __init__.py

__author__ = ', '.join([
    'Michael E. Vinyard'
])
__email__ = ', '.join([
    'vinyard@g.harvard.edu',
])


from . import _plotting  as pl
from . import _utilities as ut
from . import _tools as tl
from . import _single_cell as sc
from . import _testing_and_development as dev
from . import _data

# matplotlib defaults
import matplotlib
import matplotlib.font_manager
font = {"size": 12}
matplotlib.rc('font', **font)
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"
