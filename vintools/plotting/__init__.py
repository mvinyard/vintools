# plotting

from .fig_presets import single_fig_presets as presets
from .make_data_gif import make_gif as make_gif


from .color_palettes import vin_colors
from .color_palettes import shareseq_colors as share_seq

from .share_seq_df_generation import shareseq_palette_df

# from ._histogram import _plot_histogram as histogram
from ._delete_spines import _delete_spines as delete_spines
from ._histogram_2d_component_plot import _histogram_2d_component_plot as hist2d_component_plot