import os
import pandas as pd

def _get_pypi_package_loc(package='vintools', tmp_info_file="info.txt"):
    
    """"""

    os.system('pip show {} > {}'.format(package, tmp_info_file))
    df = pd.read_table("info.txt", sep = ": ", engine='python')
    package_loc = df.loc[df.Name == 'Location'].vintools.values.astype(str)[0]
    os.remove(tmp_info_file)
    
    return package_loc