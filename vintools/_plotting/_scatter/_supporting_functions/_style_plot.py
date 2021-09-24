
from ..._modify_ax_spines import _modify_ax_spines as ax_spines

def _spine_presets_single_ax(ax):

    """"""

    spines = ax_spines(ax)
    spines.delete(["top", "right"])
    spines.set_color("dimgrey")
    spines.set_position(position_type="outward", amount=20)
    
def _spine_presets(AxesDict):

    for row in AxesDict.keys():
        for ax in AxesDict[row].values():
            _spine_presets_single_ax(ax)
            
def _add_grid(AxesDict):
    
    """Add grids to each subplot"""
    
    for row in AxesDict.keys():
        for ax in AxesDict[row].values():
            ax.grid(c="grey", alpha=0.5, zorder=0)            
            