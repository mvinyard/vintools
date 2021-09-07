
# package imports #
# --------------- #
import numpy as np


def _get_data_bounds(data_2d, dim_1="x", dim_2="y", round_decimal=1):

    """

    Notes:
    ------
    data_2d should be an embedding or 2d data.
    """

    DataBounds = {}
    DataBounds[dim_1] = {}
    DataBounds[dim_2] = {}

    DataBounds[dim_1]["min"], DataBounds[dim_1]["max"] = (
        np.round(data_2d[:, 0].min(), decimals=round_decimal),
        np.round(data_2d[:, 0].max(), decimals=round_decimal),
    )
    DataBounds[dim_2]["min"], DataBounds[dim_2]["max"] = (
        np.round(data_2d[:, 1].min(), decimals=round_decimal),
        np.round(data_2d[:, 1].max(), decimals=round_decimal),
    )

    return DataBounds