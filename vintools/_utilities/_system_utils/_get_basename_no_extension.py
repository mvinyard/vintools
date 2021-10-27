
import os

def _get_basename_no_extension(path):

    try:
        return os.path.basename(path).split(".")[0]
    except:
        return os.path.basename(path)
