"""
AverageSpectrum(timeSeries, fsamp, numWindows, window) produces a properly
normalized power spectrum , P(f), of timeSeries, sampled at fsamp, using the
window function given by the user, and averaged over numWindows of windows.

Jeffrey Moffitt
November 15, 2004

Modified to allow smoothing of data: simple smoothing of outlier data
points
081215, mjc
"""
import numpy as np
from OneSidedSpectrum import OneSidedSpectrum
from nowindow import nowindow


def AverageSpectrum(timeSeries, fsamp, numWindows=1, window=nowindow):
    timeSeries = np.asarray(timeSeries, dtype=float).ravel()
    numWindows = int(numWindows)
    seglen = len(timeSeries) // numWindows

    index1 = 1
    index2 = seglen
    sumSpectrum, f = OneSidedSpectrum(timeSeries[index1 - 1:index2], fsamp, window)
    c = sumSpectrum.size
    tempSpectrum = np.zeros((numWindows, c))
    tempSpectrum[0, :] = sumSpectrum
    for i in range(2, numWindows + 1):
        index1 = seglen * (i - 1) + 1
        index2 = seglen * i
        tempSpectrum[i - 1, :], _ = OneSidedSpectrum(timeSeries[index1 - 1:index2], fsamp, window)
    sumSpectrum = np.sum(tempSpectrum, axis=0)
    spectrum = sumSpectrum / numWindows

    spectrumraw = spectrum
    return spectrum, spectrumraw, f
