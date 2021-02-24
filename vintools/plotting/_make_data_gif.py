from ._fig_presets import _single_fig_presets as presets
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import Image as nb_display
import glob
import os


def img_dir(dir_name):

    """Makes image directory if needed."""

    if os.path.isdir(dir_name) == False:
        os.mkdir(dir_name)


def plot_fig_single_pt(data, save_dir="./imgs/"):

    """Data must be an (N x 2) numpy array"""

    img_dir(save_dir)

    for i, sample in enumerate(data):

        presets("EMT PCA: One Trajectory", x_lab="PC-1", y_lab="PC-2")
        savename = save_dir + str(str(i).zfill(4)) + ".png"
        plt.xlim(-0.25, 0.75)
        plt.ylim(-0.035, 0.065)
        if i > 0:
            plt.scatter(
                data[: i - 1, 0], data[: i - 1, 1], c="lightgrey", alpha=0.75,
            )
        plt.scatter(data[i, 0], data[i, 1], c="purple")
        plt.savefig(savename)
        plt.close()


def _make_gif(data, save_dir="./imgs/", show=True, gif_savename="gif.gif"):

    """Data must be an (N x 2) numpy array"""
    
    from PIL import Image

    plot_fig_single_pt(data, save_dir)

    frames = []
    imgfiles = save_dir + "*.png"
    imgs = np.sort(glob.glob(imgfiles))
    for i in imgs:
        new_frame = Image.open(i)
        frames.append(new_frame)

    frames[0].save(
        gif_savename,
        format="GIF",
        append_images=frames[1:],
        save_all=True,
        duration=1,
        loop=0,
    )

    if show == True:
        with open(gif_savename, "rb") as f:
            display(nb_display(data=f.read(), format="png"))