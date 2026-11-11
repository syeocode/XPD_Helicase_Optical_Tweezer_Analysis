"""
OneSidedSpectrum(timeSeries, fsamp, window) produces a properly
normalized power spectrum , P(f), of timeSeries, sampled at fsamp, using the
window function given by the user.

Jeffrey Moffitt
November 15, 2004
"""
import numpy as np


def OneSidedSpectrum(timeSeries, fsamp, window):
    timeSeries = np.asarray(timeSeries, dtype=float).ravel()
    timeSeries = timeSeries - np.mean(timeSeries)

    # New Code 6/29/06  Make timeSeries an even number of points
    timeSeries = timeSeries[:(len(timeSeries) // 2) * 2]

    # New Code 8/3/06
    N = len(timeSeries)

    w = np.asarray(window(N), dtype=float).ravel()

    windowedSeries = w * timeSeries
    windowSum = N * np.sum(w * w)
    transform = np.fft.fft(windowedSeries)

    power = np.real(transform * np.conj(transform))

    half = N // 2
    oneSidedPower = np.zeros(half + 1)

    # New Code 8/3/06
    oneSidedPower[0] = power[0]                        # Zero frequency point
    oneSidedPower[1:half] = 2 * power[1:half]          # Positive and negative frequencies
    oneSidedPower[half] = 2 * power[half]              # fNyquist

    df = fsamp / N
    f = np.arange(half + 1) * df  # corrected 8/3/06 YC.

    Spectrum = oneSidedPower / (df * windowSum)  # Integrating the power spectrum will return
                                                 # the observed RMS.
    return Spectrum, f
