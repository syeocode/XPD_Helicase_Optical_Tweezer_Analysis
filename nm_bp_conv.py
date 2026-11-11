import numpy as np
import shared_globals as G
from XWLCContour import XWLCContour


def nm_bp_conv(ext_nm, F, fkhp):
    """fkhp 'fork' or 'hairpin'.  Uses the global WLC_param (shared_globals)."""
    WLC_param = G.WLC_param
    F = np.atleast_1d(np.asarray(F, dtype=float))

    bp = []
    nt = []
    fconv = []

    for ii in range(len(F)):
        bp.append(1 * WLC_param.hds * XWLCContour(F[ii], WLC_param.Pds, WLC_param.Sds, WLC_param.kT))  # nm/bp
        nt.append(1 * WLC_param.hss * XWLCContour(F[ii], WLC_param.Pss, WLC_param.Sss, WLC_param.kT))  # nm/nt
    bp = np.array(bp, dtype=float)
    nt = np.array(nt, dtype=float)

    if np.isnan(np.mean(bp)):
        # get indices of NaN's
        indx_bp_NaN = np.isnan(bp)
        indx_nt_NaN = np.isnan(nt)
        # get mean values without NaN's
        mean_bp = np.mean(bp[~indx_bp_NaN])
        mean_nt = np.mean(nt[~indx_nt_NaN])
        # replace NaN's with mean values
        bp[indx_bp_NaN] = mean_bp
        nt[indx_nt_NaN] = mean_nt

    # NOTE: as in the original MATLAB code this is outside the loop, so only
    # the last force value (ii = end) is used for the conversion factor.
    ii = len(F) - 1
    # nm-bp conversion factor for the fork construct
    if fkhp == 'fork':
        fconv.append(nt[ii] - bp[ii])
    # nm-bp conversion factor for the hairpin construct
    elif fkhp == 'hairpin':
        fconv.append(2 * nt[ii])
    fconv = np.array(fconv)

    ext_bp = np.asarray(ext_nm, dtype=float) / (fconv[0] if fconv.size == 1 else fconv)
    return ext_bp
