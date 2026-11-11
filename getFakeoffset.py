import numpy as np
import shared_globals as G
from ConstructConstants import ConstructConstants
from XWLCContour import XWLCContour
from matlab_helpers import fminsearch


def getFakeoffset(DNAext, FM, fitrange, fake_offset_0, construct=None):
    """
    190117 Steve
    Get fakeoffset for FEC by :
    1) Get force array of data within fitrange
    2) Get ext from data corresponding to 1)
    3) Get ext from dmodel corresponding to 1)
    4) Minimize the difference.

    Uses the global WLC_param. In MATLAB `construct` was also declared global;
    here the argument is used, falling back to shared_globals.construct.
    """
    if construct is None:
        construct = G.construct
    WLC_param = G.WLC_param
    ConstructConst = ConstructConstants(construct)

    DNAext = np.asarray(DNAext, dtype=float)
    FM = np.asarray(FM, dtype=float)

    fkoffsetidx = np.where((FM > fitrange[0]) & (FM < fitrange[1]))[0]  # get FM index within the fitrange.
    fkoffset_Farr = FM[fkoffsetidx]  # get FM force array that's within the fitrange
    DNAext_offset = DNAext[fkoffsetidx]  # get DNAext corresponding to the fkoffsetidx

    # get WLC extension using fkoffset_Farr ( force array of the data. )
    bp = 1 * WLC_param.hds * XWLCContour(fkoffset_Farr, WLC_param.Pds, WLC_param.Sds, WLC_param.kT)  # nm/bp
    nt = 1 * WLC_param.hss * XWLCContour(fkoffset_Farr, WLC_param.Pss, WLC_param.Sss, WLC_param.kT)  # nm/nt
    # closed model
    WLC_ext_fit = ConstructConst.ssDNA1 * nt + ConstructConst.dsDNA1 * bp

    # find o where the difference is minized.
    fakeoffset = fminsearch(lambda x: np.sum((WLC_ext_fit - (DNAext_offset - x))**2), fake_offset_0)
    return fakeoffset
