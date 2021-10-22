from ._supporting_functions._construct_plot import _construct_plot
from ._supporting_functions._style_plot import _spine_presets, _add_grid
from ._supporting_functions._plot_data import _plot_data, _plot_simple

import matplotlib

font = {"size": 12}
matplotlib.rc(font)
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"

import matplotlib.pyplot as plt


class ScatterPlot:
    def __init__(self, tight=True, grid=True):

        """Instantiates a scatterplot"""

    def construct_layout(
        self,
        nplots,
        ncols=4,
        figsize_width=1,
        figsize_height=1,
        figsize=False,
        grid_hspace=0,
        width_ratios=False,
    ):

        """
        Setup figure layout.

        nplots

        ncols

        figsize_width

        figsize_height
        """

        self.nplots = nplots
        self.ncols = ncols

        if figsize:
            self.figsize_width = self.figsize_height = figsize
        else:
            self.figsize_width = figsize_width
            self.figsize_height = figsize_height

        self.grid_hspace = grid_hspace
        self.width_ratios = width_ratios

        self.fig, self.AxesDict = _construct_plot(
            nplots=self.nplots,
            ncols=self.ncols,
            figsize_width=self.figsize_width,
            figsize_height=self.figsize_height,
            grid_hspace=self.grid_hspace,
            width_ratios=self.width_ratios,
        )

    def style(self, grid=True):

        _spine_presets(self.AxesDict)
        if grid:
            _add_grid(self.AxesDict)

    def plot_data(self, adata, embedding, variables_to_plot):

        """ """
        self.plots = _plot_data(self.AxesDict, adata, embedding, variables_to_plot)

    def simple(self, x, y, title=None, color="grey"):

        self.plots = _plot_simple(self.AxesDict, x, y, title, color)


def _scatter(
    adata=False,
    x=False,
    y=False,
    title=None,
    color="grey",
    variables_to_plot=None,
    embedding="X_SPRING",
    show=True,
):

    """
    adata
       
    variables_to_plot
    
    embedding

    """

    if variables_to_plot == None:
        nplots = 1
    else:
        nplots = len(variables_to_plot)

    sc = ScatterPlot()
    sc.construct_layout(nplots=nplots)
    sc.style()

    if adata:
        plots = sc.plot_data(adata, embedding, variables_to_plot)
    else:
        plots = sc.simple(x, y, title, color)
    if show:
        plt.show()
