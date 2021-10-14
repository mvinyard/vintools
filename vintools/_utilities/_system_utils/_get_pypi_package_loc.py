import os
import string
import numpy as np
import pandas as pd


def _get_pypi_package_loc(package="vintools"):

    """"""

    signature = "".join(np.random.choice(list(string.ascii_lowercase), 6))
    tmp_info_file = "info_{}.txt".format(signature)

    os.system("pip show {} > {}".format(package, tmp_info_file))
    df = pd.read_table(tmp_info_file, sep=": ", engine="python")
    package_loc = df.loc[df.Name == "Location"].vintools.values.astype(str)[0]
    os.remove(tmp_info_file)

    return package_loc
