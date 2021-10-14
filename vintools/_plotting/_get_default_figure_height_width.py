import matplotlib


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
