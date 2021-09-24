
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

def _calculate_nrows(nplots, ncols):

    return math.ceil(nplots / ncols)

def _get_default_figure_height_width():

    """
    Return default height and width of matplotlib figures.

    Parameters:
    -----------
    None

    Returns:
    --------
    DefaultFigsizeDict
        two key-value pairs of height and width.
        type: Dict
    """

    DefaultFigsizeDict = {}

    default_figsize_mpl = matplotlib.rcParams["figure.figsize"]
    DefaultFigsizeDict["width"], DefaultFigsizeDict["height"] = (
        default_figsize_mpl[0],
        default_figsize_mpl[1],
    )
    
    return DefaultFigsizeDict

def _setup_fig(ncols, nrows, figsize_width, figsize_height):

    """

    Parameters:
    -----------
    ncols
        Number of columns in the figure.
        type: int

    nrows
        Number of rows in the figure.
        type: int

    figsize_width
        Scaler adjustment of figure width
        default: 1
        type: float

    figsize_height
        Scaler adjustment of figure height
        default: 1
        type: float

    Returns:
    --------
    fig
        type: matplotlib.figure.Figure

    Notes:
    ------
    """

    default_figsize = _get_default_figure_height_width()
    h_def = default_figsize["height"]
    w_def = default_figsize["width"]

    fig_width = w_def * ncols * figsize_width
    fig_height = h_def * nrows * figsize_height

    fig = plt.figure(figsize=(fig_width, fig_height))

    return fig

def _construct_plot(nplots, ncols=4, figsize_width=1, figsize_height=1):

    """
    Creates Axes for each desired plot.

    Parameters:
    -----------
    nplots

    ncols
        Number of columns. 
        default: 4
        type: int

    Returns:
    --------

    Notes:
    ------

    """
    
    nrows = _calculate_nrows(nplots, ncols)
    fig = _setup_fig(ncols, nrows, figsize_width, figsize_height)
    gridspec = GridSpec(nrows, ncols, width_ratios=np.ones(ncols))
    
    plot_count = 0
    AxesDict = {}

    for ax_i in range(nrows):
        AxesDict[ax_i] = {}
        for ax_j in range(ncols):
            plot_count += 1
            AxesDict[ax_i][ax_j] = fig.add_subplot(gridspec[ax_i, ax_j])
            if plot_count >= nplots:
                break

    return fig, AxesDict