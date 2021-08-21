
# package imports #
# --------------- #
# import matplotlib
# import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar

# local imports #
# ------------- #
from ._modify_ax_spines import _modify_ax_spines
from ._set_minimal_ticks import _set_minimal_ticks

# # matplotlib presets #
# # ------------------ #
# font = {"size": 12}
# matplotlib.rc(font)
# matplotlib.rcParams["font.sans-serif"] = "Arial"
# matplotlib.rcParams["font.family"] = "sans-serif"

def _histogram_2d_component_plot(
    data,
    n_bins=20,
    cmap="Purples",
    histcolor="#63439C",
    histalpha=1,
    figsize=(7, 6),
    title_fontsize=16, 
    label_fontsize=14,
    suptitle=False,
    save_path=False,
):

    """"""

    x_component, y_component = data[:, 0], data[:, 1]

    # construct plot surface
    fig = plt.figure(figsize=figsize)
    gridspec = GridSpec(
        nrows=2,
        ncols=3,
        height_ratios=[0.15, 1],
        width_ratios=[0.2, 1, 0.075],
        wspace=0.2,
        hspace=0.1,
    )

    x_component_histogram = fig.add_subplot(gridspec[0, 1])
    x_component_bins = x_component_histogram.hist(
        x_component,
        bins=n_bins,
        color=histcolor,
        alpha=histalpha,
    )
    x_component_histogram.grid(True)
#     x_component_histogram.set_xticks([])

    y_component_histogram = fig.add_subplot(gridspec[1, 0])
    y_component_bins = y_component_histogram.hist(
        y_component,
        orientation="horizontal",
        bins=n_bins,
        color=histcolor,
        alpha=histalpha,
    )
    y_component_histogram.grid(True)

    y_component_histogram.set_xlim(y_component_histogram.get_xlim()[::-1])

    histogram_2d = fig.add_subplot(gridspec[1, 1])
    bins_2d = histogram_2d.hist2d(x_component, y_component, bins=n_bins, cmap=cmap)

    #     histogram_2d.set_yticks([bins_2d[2].min() + 1, bins_2d[2].max() - 1])

    x_tick_min = int(bins_2d[1].min())
    x_tick_max = int(bins_2d[1].max())

    histogram_2d.set_xticks([x_tick_min, x_tick_max])
    #

    y_tick_min = int(bins_2d[2].min())
    y_tick_max = int(bins_2d[2].max())

    histogram_2d.set_yticks([y_tick_min, y_tick_max])
    histogram_2d.yaxis.tick_right()
#     y_component_histogram.set_yticks([])

    hist2d_spines = _modify_ax_spines(histogram_2d)
    hist2d_spines.set_color("grey")
#     spines.delete(select_spines=["top", "right"])
    hist2d_spines.set_position(position_type="axes", amount=-0.05)
    _set_minimal_ticks(histogram_2d, x_component, y_component)
    
#     for ax in [x_component_histogram, y_component_histogram, histogram_2d]:
#         spines = _modify_ax_spines(ax)
#         spines.delete()
        
    
    

    hist2d_im = histogram_2d.imshow(data, cmap=cmap)
    colorbar_ax = fig.add_subplot(gridspec[1:, 2])

    cb = Colorbar(ax=colorbar_ax, mappable=hist2d_im, ticklocation="right")
    cb.outline.set_visible(False)
    if suptitle:
        plt.suptitle(suptitle, fontsize=title_fontsize, x=.55)
    plt.tight_layout()
    if save_path:
        out = plt.savefig(save_path)
    plt.show()
