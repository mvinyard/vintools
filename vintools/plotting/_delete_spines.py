def _delete_spines(ax, spines=["top", "right", "bottom", "left"]):

    "https://newbedev.com/how-to-remove-frame-from-matplotlib-pyplot-figure-vs-matplotlib-figure-frameon-false-problematic-in-matplotlib"

    for spine in spines:
        ax.spines[spine].set_visible(False)