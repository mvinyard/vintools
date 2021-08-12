
import numpy as np

def _set_minimal_ticks(ax, x, y, round_decimal=2):

    """"""

    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    x_mid, y_mid = np.mean([x.min(), x.max()]), np.mean([y.min(), y.max()])
    x_mid_low, x_mid_high = (
        np.mean([x_mid, x.min()]),
        np.mean([x_mid, x.max()]),
    )
    y_mid_low, y_mid_high = (
        np.mean([y_mid, y.min()]),
        np.mean([y_mid, y.max()]),
    )
    x_mid, y_mid = x_mid, y_mid

    x_ticks = np.round(np.array([x_min, x_mid_low, x_mid, x_mid_high, x_max]), decimals=round_decimal)
    y_ticks = np.round(np.array([y_min, y_mid_low, y_mid, y_mid_high, y_max]), decimals=round_decimal)

    ax.get_xaxis().set_ticks(x_ticks)
    ax.get_yaxis().set_ticks(y_ticks)