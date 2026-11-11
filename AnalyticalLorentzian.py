"""
AnalyticalLorentzian
[fc, D, sfc, sD, chi] = AnalyticalLorentzian(spec, df, Tmsr, fstart)
This function is based on  K. Berg-Sorensen, H. Flyvberg, Rev. Sci.
Inst., 75, 594 (2004).  It fits a lorentzian to a positive frequency
power spectra, returning the diffusion constant and the corner frequency
and the uncertainty in each of these.

Jeffrey Moffitt
Corrected chi 2/4/06
"""
import numpy as np


def AnalyticalLorentzian(spec, df, Tmsr, fstart=0):
    # MATLAB returns NaN/Inf silently (e.g. for a flat Sum signal); do the same
    with np.errstate(divide='ignore', invalid='ignore'):
        return _AnalyticalLorentzian(spec, df, Tmsr, fstart)


def _AnalyticalLorentzian(spec, df, Tmsr, fstart=0):
    spec = np.asarray(spec, dtype=float).ravel()

    # Compute useful quantities
    f = np.arange(len(spec)) * df + fstart

    S00 = len(spec)
    S01 = np.sum(spec)
    S02 = np.sum(spec * spec)
    S11 = np.sum(spec * f * f)
    S12 = np.sum(spec * spec * f * f)
    S22 = np.sum(spec * spec * f * f * f * f)

    fc = np.sqrt(complex((S01 * S22 - S11 * S12) / (S11 * S02 - S01 * S12)))
    fc = fc.real if fc.imag == 0 else fc
    D = np.pi**2 * (S02 * S22 - S12**2) / (S11 * S02 - S01 * S12)

    chi = S00 - (S01**2 * S22 + S11**2 * S02 - 2 * S01 * S11 * S12) / (S02 * S22 - S12**2)

    x1 = fstart / fc
    x2 = (fstart + len(spec) * df) / fc

    u = 2 * x2 / (1 + x2**2) - 2 * x1 / (1 + x1**2) + 2 * np.arctan((x2 - x1) / (1 + x2 * x1))

    v = 4 / (x2 - x1) * (np.arctan((x2 - x1) / (1 + x2 * x1)))**2

    sfc = np.sqrt(fc / ((u - v) * Tmsr))

    sD = D * np.sqrt(u / (fc * Tmsr * (u - v) * (x2 - x1)))

    return fc, D, sfc, sD, chi
