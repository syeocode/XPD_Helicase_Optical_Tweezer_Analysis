import glob
import os
import numpy as np
from scipy.io import loadmat


def getBetaforOT(Pathtofolder):
    """
    Get beta from the Folder. The beta will be stored in filename fixedxyres_YYMMDD.mat
    231101  SY
    """
    f = sorted(glob.glob(os.path.join(Pathtofolder, 'fixedxyres*.mat')))
    if len(f) == 0:
        raise FileNotFoundError(f'No fixedxyres*.mat file in {Pathtofolder}')
    beta = loadmat(f[0])
    beta = np.asarray(beta['beta'], dtype=float)[0:2, 0:2]
    return beta
