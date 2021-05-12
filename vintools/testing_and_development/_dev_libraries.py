def _load_development_libraries():
    
    """
    Assigns global variables to packages used in the development of nodescape such that developing with common packages can be done cleanly.
    
    """
    
    global glob
    global torch
    global time
    global np
    global pd
    global plt
    global nn
    global a
    global os
    global optim
    global odeint
    global sp
    global PCA
    global v
    
    import glob
    import torch as torch
    from torchdiffeq import odeint
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import torch.nn as nn
    import torch.optim as optim
    import anndata as a
    import os    
    import time
    import scipy as sp
    from sklearn.decomposition import PCA
    
    return odeint, torch, np, pd, plt, nn, a, os, glob, time, optim, sp, PCA