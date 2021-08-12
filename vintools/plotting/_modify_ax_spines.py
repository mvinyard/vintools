def _choose_spines(self, select_spines):
    if not select_spines:
        spines = self.spines
    else:
        spines = select_spines

    return spines


class _modify_ax_spines:
    def __init__(self, ax):

        self.ax = ax
        self.spines = ["left", "right", "bottom", "top"]

    def set_color(self, color, select_spines=False):

        """
        Set the position of each spine.

        Parameters:
        -----------
        color
            type: str

        select_spines
            default: False
            type: str (or default, bool: False)

        Returns:
        --------
        None
            modifies ax

        Notes:
        ------
        """

        spines = _choose_spines(self, select_spines)
        for spine in spines:
            self.ax.spines[spine].set_color(color)

    def delete(self, select_spines=False):

        """
        Parameters:
        -----------
        select_spines
            default: False
            type: str (or default, bool: False)

        Returns:
        --------
        None
            modifies ax

        Notes:
        ------
        (1) source: https://newbedev.com/how-to-remove-frame-from-matplotlib-pyplot-figure-vs-matplotlib-figure-frameon-false-problematic-in-matplotlib
        """

        spines = _choose_spines(self, select_spines)
        for spine in spines:
            self.ax.spines[spine].set_visible(False)

    def set_position(self, position_type, amount, select_spines=False):

        """
        Set the position of each spine.

        Parameters:
        -----------
        position_type
            type:str

        amount
            type: float

        Returns:
        --------
        None
            modifies ax

        Notes:
        ------
        (1) From the matplotlib documentation:
            https://matplotlib.org/stable/api/spines_api.html#matplotlib.spines.Spine.set_position

            Spine position is specified by a 2 tuple of (position type, position type). The position types are:
            'outward': place the spine out from the data area by the specified number of points. (Negative values place the spine inwards.)
            'axes': place the spine at the specified Axes coordinate (0 to 1).
            'data': place the spine at the specified data coordinate.
        (2) Additionally, shorthand notations define a special positions:

            'center' -> ('axes', 0.5)
            'zero' -> ('data', 0.0)
        """

        spines = _choose_spines(self, select_spines)
        for spine in spines:
            self.ax.spines[spine].set_position((position_type, amount))