def single_fig_presets(title, x_lab, y_lab, size=(10,8), mute_x_axis=False):

    """
    presets for one single figure to look nice
    """

    import matplotlib.pyplot as plt

    fig = plt.figure(figsize = size)
    ax = fig.add_subplot(1,1,1)
    ax.set_title(title, fontsize = 20)
    ax.set_xlabel(x_lab, fontsize = 15)
    ax.set_ylabel(y_lab, fontsize = 15)
    ax.spines['left'].set_linewidth(3)
    ax.spines['bottom'].set_linewidth(3)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    if mute_x_axis == True:
        ax.axes.get_xaxis().set_visible(False)
        
    return fig, ax