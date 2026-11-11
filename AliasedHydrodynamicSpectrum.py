import numpy as np
from HydrodynamicSpectrum import HydrodynamicSpectrum


def AliasedHydrodynamicSpectrum(f, para, fsamp, n):
    """
    spec = AliasedHydrodynamicSpectrum(f, [fc, D, fv, fm], fsamp, n)

    spec is the sum of a Hydrodynamic Spectrum evaluated at f + i*fsamp where i=-n:n
    This produces a spectrum that well approximates the effects of the
    inherent aliasing produced by sampling at fsamp.  Note that no
    anti-aliasing filtering is assumed.
    Based on K. Berg-Sorensen, H. Flyvbjerg. Rev Sci Inst, 75 (2004)

    Jeffrey Moffitt
    January 24, 2006
    jmoffitt@berkeley.edu
    """
    f = np.asarray(f, dtype=float)
    spec = 0 * f  # Kludge so that spec has the same shape as f
    for i in range(-n, n + 1):
        spec = spec + HydrodynamicSpectrum(f + i * fsamp, para)
    return spec
