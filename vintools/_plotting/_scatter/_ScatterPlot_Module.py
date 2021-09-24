
from ._supporting_functions._construct_plot import _construct_plot
from ._supporting_functions._style_plot import _spine_presets, _add_grid
from ._supporting_functions._plot_data import _plot_data

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
        self, nplots, ncols=4, figsize_width=1, figsize_height=1, figsize=False
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

        self.fig, self.AxesDict = _construct_plot(
            nplots=self.nplots, ncols=self.ncols, figsize_width=self.figsize_width, figsize_height=self.figsize_height
        )
        
    def style(self):
        
        _spine_presets(self.AxesDict)
        _add_grid(self.AxesDict)
        
    def plot_data(self, adata, embedding, variables_to_plot):
        
        """ """        
        _plot_data(self.AxesDict, adata, embedding, variables_to_plot)
        
        
def _scatter(adata, variables_to_plot, embedding="X_SPRING"):
    
    """
    adata
       
    variables_to_plot
    
    embedding

    """

    sc = ScatterPlot()
    sc.construct_layout(nplots=len(variables_to_plot))
    sc.style()
    sc.plot_data(adata, embedding, variables_to_plot)
    plt.show()