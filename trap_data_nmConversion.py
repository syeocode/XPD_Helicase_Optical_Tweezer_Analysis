"""
Originally by Monika; updated by AVT

230802: adding TrapSep_X and TrapSep_Y as variables
230818: saving above as fields
230910: implementing beta and gamma conversion for OT.
-- beta: 2*2 matrix of V-to-pixel conversion factors, based on bead grid images
-- gamma: nm/pixel conversion factor from image of 10 um grid
"""
import copy
import numpy as np


def trap_data_nmConversion(trapData, param, beta=None):
    trapData = copy.copy(trapData)  # MATLAB passes structs by value
    if beta is None:
        beta = np.eye(2)  # dummy identity matrix

    flconv = 123
    otconvX = 883.59  # nm/V updated on 06-29-2018
    otconvY = 508.93  # nm/V updated on 06-29-2018

    if param.instrument == 0:
        trapData.TrapSep = (np.asarray(trapData.trappos2) - np.asarray(trapData.trappos1)) * flconv
    elif param.instrument == 1:
        trapData.TrapSep = np.sqrt(((np.asarray(trapData.trappos1) - trapData.t2x) * otconvX)**2 +
                                   ((np.asarray(trapData.trappos2) - trapData.t2y) * otconvY)**2)
    return trapData
