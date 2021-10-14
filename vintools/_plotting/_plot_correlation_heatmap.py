
from ._plot_utils._get_default_matplotlib_figure_width_height import _get_default_matplotlib_figure_width_height

from matplotlib.gridspec import GridSpec
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import os

os.system("rm ~/.cache/matplotlib -rf")

font = {"size": 12}
matplotlib.rc(font)
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"


def _plot_correlation_heatmap(
    df, title="Sample Correlation", figsize=1, title_y=1.15, title_x=0, savename=False
):

    DefaultFigsizeDict = _get_default_figure_height_width()
    default_h, default_w = DefaultFigsizeDict["height"], DefaultFigsizeDict["width"]

    fig = plt.figure(figsize=(default_w * figsize, default_h * figsize))
    gridspec = GridSpec(1, 1)
    ax = fig.add_subplot(gridspec[0, 0])
    im = ax.imshow(df.values.astype(float), cmap="plasma", vmin=0, vmax=1)
    plt.xticks(
        range(len(df.columns.values.astype(str))),
        df.columns.values.astype(str),
        rotation=60,
        ha="left",
    )
    ax.xaxis.tick_top()
    plt.yticks(
        range(len(df.columns.values.astype(str))),
        df.columns.values.astype(str),
        rotation=0,
    )
    plt.colorbar(im, shrink=0.6, aspect=30)
    plt.title(title, y=title_y, x=title_x, fontsize=16)

    if savename:
        figsavename = savename + ".png"
        plt.savefig(figsavename, bbox_inches="tight")
    plt.tight_layout()
    plt.show()

    return figsavename
