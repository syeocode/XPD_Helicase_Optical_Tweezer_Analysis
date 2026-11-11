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


def trap_data_nmConversion_newbeta(trapData, param, beta=None):
    trapData = copy.copy(trapData)  # MATLAB passes structs by value
    if beta is None:
        beta = np.eye(2)  # dummy identity matrix
    beta = np.asarray(beta, dtype=float)

    flconv = 123
    otconvX = 883.59  # updated on 06-29-2018
    otconvY = 508.93  # updated on 06-29-2018
    ot_gamma = 144.8  # nm/pixel, average for x and y (based on 2018 grid image analysis by SY)

    if param.instrument == 0:
        # finish when gamma value for fleezers available
        trapData.TrapSep = (np.asarray(trapData.trappos2) - np.asarray(trapData.trappos1)) * flconv
        trapData.TrapSep_X = trapData.TrapSep
        trapData.TrapSep_Y = trapData.TrapSep * 0
    elif param.instrument == 1:
        if not np.array_equal(beta, np.eye(2)):
            gamma_X = 147.2863  # ot_gamma; %148.0282 SY code
            gamma_Y = 152.5629  # -ot_gamma;  % pixel y-axis is flipped; flip back (?) %152.2694 SY code
        else:
            gamma_X = otconvX
            gamma_Y = otconvY
        Trap_Sep_pixels = beta @ np.vstack([np.asarray(trapData.trappos1).ravel() - trapData.t2x,
                                            np.asarray(trapData.trappos2).ravel() - trapData.t2y])
        trapData.TrapSep_X = Trap_Sep_pixels[0, :] * gamma_X
        trapData.TrapSep_Y = Trap_Sep_pixels[1, :] * gamma_Y
        trapData.TrapSep = np.sqrt(trapData.TrapSep_X**2 + trapData.TrapSep_Y**2)

    return trapData
