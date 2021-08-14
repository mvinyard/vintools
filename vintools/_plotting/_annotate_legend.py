
# package imports #
# --------------- #
import matplotlib
import matplotlib.font_manager
import matplotlib.pyplot as plt

# matplotlib presets #
# ------------------ #
font = {"size": 12}
matplotlib.rc(font)
matplotlib.rcParams["font.sans-serif"] = "Arial"
matplotlib.rcParams["font.family"] = "sans-serif"


def _annotate_legend(ax, 
    markerscale=3, fontsize=12, loc=2, **kwargs
):

    """"""
    
    plt.legend(
        markerscale=markerscale,
        edgecolor="w",
        fontsize=fontsize,
        handletextpad=None,
        loc=loc,
    )