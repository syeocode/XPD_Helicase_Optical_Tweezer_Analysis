import numpy as np


def HydrodynamicSpectrum(f, para):
    """spec = HydrodynamicSpectrum(f, [fc, D, fv, fm])"""
    fc = para[0]
    D = para[1]
    fv = para[2]
    fm = para[3]
    f = np.abs(np.asarray(f, dtype=float))  # Kludge!
    spec = D / np.pi**2 * (1 + np.sqrt(f / fv)) / (
        (fc - f * np.sqrt(f / fv) - (f**2) / fm)**2 + (f + f * np.sqrt(f / fv))**2)
    return spec
