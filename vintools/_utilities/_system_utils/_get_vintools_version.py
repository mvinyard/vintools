def _get_vintools_version(setup_path="/home/mvinyard/software/vintoos/setup.py"):
    
    """reports version without circular import of package"""
    
    with open(setup_path, 'r') as setup_py:
        for n, line in enumerate(setup_py):
            l = line.rstrip()
            if l.startswith('    version'):
                version = line.split('"')[-2]
                return version