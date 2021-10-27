
# _fix_path_for_glob.py
__module_name__ = "_fix_path_for_glob.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])


def _fix_path_for_glob(path):

    if path.endswith("/*"):
        pass
    elif path.endswith("/"):
        path = path + "*"
    else:
        path = path + "/*"

    return path
