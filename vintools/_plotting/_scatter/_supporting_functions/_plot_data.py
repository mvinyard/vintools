
import numpy as np
from ...color_palettes import vin_colors

def _plot_categorical(ax, x, y, df, variable):

    """"""

    colors = vin_colors()

    df = df.reset_index()

    for i, label in enumerate(np.sort(df[variable].unique())):

        label_idx = df.loc[df[variable] == label].index.astype(int)
        x_, y_ = x[label_idx], y[label_idx]
        ax.scatter(x_, y_, c=colors[i], zorder=1000, label=label)


def _get_GEX(adata, gene):

    GEX_values = adata[:, gene].X.toarray()

    return GEX_values


def _plot_continuous(ax, adata, x, y, variable):

    if variable in adata.var_names:
        color_values = _get_GEX(adata, variable)
    else:
        color_values = adata.obs[variable].values.astype(float)

    if variable == None:
        color_values = "lightgrey"

    ax.scatter(x, y, c=color_values, zorder=1000)


def _make_subplot(ax, adata, embedding, variable=None):

    x, y = adata.obsm[embedding][:, 0], adata.obsm[embedding][:, 1]

    if (type(variable) is str) and (not variable in adata.var_names):
        _plot_categorical(ax, x, y, df=adata.obs, variable=variable)
    else:
        _plot_continuous(ax, adata, x, y, variable=variable)

    ax.set_title(variable)


def _plot_data(AxesDict, adata, embedding, variables_to_plot):

    if type(variables_to_plot) == str:
        variables_to_plot = [variables_to_plot]

    for row in AxesDict.keys():
        for n_plot, ax in enumerate(AxesDict[row].values()):
            if len(variables_to_plot) == n_plot:
                break
            _make_subplot(ax, adata, embedding, variable=variables_to_plot[n_plot])
            
            
##### NON-ANNDATA ###### 

def _make_simplesubplot(ax, x, y, title, color):


    ax.scatter(x, y, c=color, zorder=1000)
    ax.set_title(title)
        
def _plot_simple(AxesDict, x, y, title, color):
    
    for row in AxesDict.keys():
        for n_plot, ax in enumerate(AxesDict[row].values()):
            _make_simplesubplot(ax, x, y, title, color)