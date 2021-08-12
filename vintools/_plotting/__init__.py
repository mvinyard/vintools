# plotting

from .fig_presets import single_fig_presets as presets
from .make_data_gif import make_gif as make_gif


from .color_palettes import vin_colors
from .color_palettes import shareseq_colors as share_seq

from ._annotate_legend import _annotate_legend as legend

from .share_seq_df_generation import shareseq_palette_df

from ._set_minimal_ticks import _set_minimal_ticks as set_minimal_ticks

# from ._histogram import _plot_histogram as histogram
from ._modify_ax_spines import _modify_ax_spines as ax_spines
from ._histogram_2d_component_plot import _histogram_2d_component_plot as hist2d_component_plot