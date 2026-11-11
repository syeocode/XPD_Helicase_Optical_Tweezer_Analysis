"""
Author: Matt Comstock
Date last modified: never (at least not by Kevin)

This function appears to have come direcly from CalibrateExact7 (or some
variant thereof). Only one bead calibration is done in this version,
however, since the two beads are handled separately in
CalibrateStrippedWrapped.

basic calibration routine that just generates the data, nothing else, for
a single set of data
090514 mjc
"""
from datetime import datetime, date
import numpy as np

from matlab_helpers import Struct, mstd
from AverageSpectrum import AverageSpectrum
from AnalyticalLorentzian import AnalyticalLorentzian
from NonLinearFit import NonLinearFit
from HydrodynamicSpectrum import HydrodynamicSpectrum
from AliasedHydrodynamicSpectrum import AliasedHydrodynamicSpectrum


def CalibrateExact_Stripped_Data_SY(param, fitmodel, data):
    fit = Struct()
    calibration = Struct()

    n = 9e-10         # [pN s/nm^2] water viscosity at 24.4C
    # n = 9e-10*1.4   # [pN s/nm^2] water viscosity at 24.4C, for 10% glycerol
    p = 1e-21         # [pN s^2/nm^4] water density
    pbead = 1.05e-21  # [pN s^2/nm^4] bead density

    beadA = param[0]
    gA = 3 * np.pi * n * beadA
    mA = np.pi / 6 * pbead * beadA**3
    fXYHigh = param[3]
    fSumHigh = param[4]
    fLow = param[5]
    avWin = int(param[6])

    # Set Hydrodynamic Parameters
    fmA = gA / (2 * np.pi * (mA + np.pi / 12 * p * beadA**3))
    fvA = n / p / (np.pi * beadA**2 / 4)

    # can determine from timing data
    fsamp = 1 / data.sampperiod

    # hardwire background subtract off
    AXSpecBkgd = 0
    AYSpecBkgd = 0
    ASumSpecBkgd = 0

    X = np.asarray(data.X, dtype=float)
    Y = np.asarray(data.Y, dtype=float)
    Sum = np.asarray(data.Sum, dtype=float)

    # Calculate Spectra
    fit.AXSpec, fit.AXSpecRaw, fit.f = AverageSpectrum(X - np.mean(X), fsamp, avWin)
    fit.AYSpec, fit.AYSpecRaw, fit.f = AverageSpectrum(Y - np.mean(Y), fsamp, avWin)
    fit.ASumSpec, fit.ASumSpecRaw, fit.f = AverageSpectrum(Sum - np.mean(Sum), fsamp, avWin)

    # Subtract background
    fit.AXSpec = fit.AXSpec - AXSpecBkgd
    fit.AYSpec = fit.AYSpec - AYSpecBkgd
    fit.ASumSpec = fit.ASumSpec - ASumSpecBkgd

    indXY = np.where((fit.f > fLow) & (fit.f < fXYHigh))[0]
    indSum = np.where((fit.f > fLow) & (fit.f < fSumHigh))[0]

    # Use Analytical Lorentizian to Guess Initial Fit Values
    # (arguments as in the original: Tmsr = fLow, fstart = 10)
    iguess_AXfc, iguess_AXD = AnalyticalLorentzian(fit.AXSpec[indXY], fit.f[1], fLow, 10)[:2]
    iguess_AYfc, iguess_AYD = AnalyticalLorentzian(fit.AYSpec[indXY], fit.f[1], fLow, 10)[:2]
    iguess_ASumfc, iguess_ASumD = AnalyticalLorentzian(fit.ASumSpec[indSum], fit.f[1], fLow, 10)[:2]

    # Fit with more sophisticated model
    model = fitmodel.lower()
    if model in ('hydrodynamicspectrum', 'hydro'):
        specFunc = lambda f, para: HydrodynamicSpectrum(f, [para[0], para[1], fvA, fmA])
        fit.AXfc, fit.AXD = NonLinearFit(fit.f[indXY], fit.AXSpec[indXY], specFunc, [iguess_AXfc, iguess_AXD])[:2]
        fit.AYfc, fit.AYD = NonLinearFit(fit.f[indXY], fit.AYSpec[indXY], specFunc, [iguess_AYfc, iguess_AYD])[:2]
        fit.ASumfc, fit.ASumD = NonLinearFit(fit.f[indSum], fit.ASumSpec[indSum], specFunc, [iguess_ASumfc, iguess_ASumD])[:2]
        predictedAX = HydrodynamicSpectrum(fit.f, [fit.AXfc, fit.AXD, fvA, fmA])
        predictedAY = HydrodynamicSpectrum(fit.f, [fit.AYfc, fit.AYD, fvA, fmA])
        predictedASum = HydrodynamicSpectrum(fit.f, [fit.ASumfc, fit.ASumD, fvA, fmA])

    elif model in ('aliasedspectrum', 'aliased'):
        # AliasedLorentzianSpectrum was not part of the converted folder
        from AliasedLorentzianSpectrum import AliasedLorentzianSpectrum
        specFunc = lambda f, para: AliasedLorentzianSpectrum(f, para, fsamp, 10)
        fit.AXfc, fit.AXD = NonLinearFit(fit.f[indXY], fit.AXSpec[indXY], specFunc, [iguess_AXfc, iguess_AXD])[:2]
        fit.AYfc, fit.AYD = NonLinearFit(fit.f[indXY], fit.AYSpec[indXY], specFunc, [iguess_AYfc, iguess_AYD])[:2]
        fit.ASumfc, fit.ASumD = NonLinearFit(fit.f[indSum], fit.ASumSpec[indSum], specFunc, [iguess_ASumfc, iguess_ASumD])[:2]
        predictedAX = AliasedLorentzianSpectrum(fit.f, [fit.AXfc, fit.AXD], fsamp, 10)
        predictedAY = AliasedLorentzianSpectrum(fit.f, [fit.AYfc, fit.AYD], fsamp, 10)
        predictedASum = AliasedLorentzianSpectrum(fit.f, [fit.ASumfc, fit.ASumD], fsamp, 10)

    elif model in ('aliasedhydrodynamicspectrum', 'aliasedhydro'):
        specFunc = lambda f, para: AliasedHydrodynamicSpectrum(f, [para[0], para[1], fvA, fmA], fsamp, 10)
        fit.AXfc, fit.AXD = NonLinearFit(fit.f[indXY], fit.AXSpec[indXY], specFunc, [iguess_AXfc, iguess_AXD])[:2]
        fit.AYfc, fit.AYD = NonLinearFit(fit.f[indXY], fit.AYSpec[indXY], specFunc, [iguess_AYfc, iguess_AYD])[:2]
        predictedAX = AliasedHydrodynamicSpectrum(fit.f, [fit.AXfc, fit.AXD, fvA, fmA], fsamp, 10)
        predictedAY = AliasedHydrodynamicSpectrum(fit.f, [fit.AYfc, fit.AYD, fvA, fmA], fsamp, 10)

    elif model in ('aliasedfilteredhydrodynamicspectrum', 'aliasedfiltered'):
        # To account for filtering of PSDs (don't use for QPD detectors!)
        # CalParameters / AliasedFilteredHydrodynamicSpectrum were not part of the converted folder
        from CalParameters import CalParameters
        from AliasedFilteredHydrodynamicSpectrum import AliasedFilteredHydrodynamicSpectrum
        calpara = CalParameters(np.mean(Sum), np.mean(Sum))
        filt = list(np.atleast_1d(calpara.AX)) + [0]
        specFunc = lambda f, para: AliasedFilteredHydrodynamicSpectrum(f, [para[0], para[1], fvA, fmA], fsamp, 10, filt)
        fit.AXfc, fit.AXD = NonLinearFit(fit.f[indXY], fit.AXSpec[indXY], specFunc, [iguess_AXfc, iguess_AXD])[:2]
        fit.AYfc, fit.AYD = NonLinearFit(fit.f[indXY], fit.AYSpec[indXY], specFunc, [iguess_AYfc, iguess_AYD])[:2]
        fit.ASumfc, fit.ASumD = NonLinearFit(fit.f[indSum], fit.ASumSpec[indSum], specFunc, [iguess_ASumfc, iguess_ASumD])[:2]
        predictedAX = AliasedFilteredHydrodynamicSpectrum(fit.f, [fit.AXfc, fit.AXD, fvA, fmA], fsamp, 10, filt)
        predictedAY = AliasedFilteredHydrodynamicSpectrum(fit.f, [fit.AYfc, fit.AYD, fvA, fmA], fsamp, 10, filt)
        predictedASum = AliasedHydrodynamicSpectrum(fit.f, [fit.ASumfc, fit.ASumD, fvA, fmA], fsamp, 10)

    else:
        # (the MATLAB code referred to undefined 'initguess_*' names here;
        #  the analytical Lorentzian guesses are what was intended)
        fit.AXfc = iguess_AXfc
        fit.AXD = iguess_AXD
        fit.AYfc = iguess_AYfc
        fit.AYD = iguess_AYD
        fit.ASumfc = iguess_ASumfc
        fit.ASumD = iguess_ASumD

        predictedAX = fit.AXD / np.pi**2 * 1 / (fit.AXfc**2 + fit.f**2)
        predictedAY = fit.AYD / np.pi**2 * 1 / (fit.AYfc**2 + fit.f**2)
        predictedASum = fit.ASumD / np.pi**2 * 1 / (fit.ASumfc**2 + fit.f**2)

    fit.predictedAX = predictedAX
    fit.predictedAY = predictedAY

    calibration.path = data.path
    calibration.file = data.file
    calibration.stamp = datetime.now()
    calibration.date = date.today()
    calibration.beadA = beadA
    calibration.fsamp = fsamp
    calibration.fXYHigh = fXYHigh
    calibration.fSumHigh = fSumHigh
    calibration.fLow = fLow

    calibration.alphaAX = np.sqrt(4.14 / (gA * fit.AXD))
    calibration.alphaAY = np.sqrt(4.14 / (gA * fit.AYD))

    calibration.kappaAX = 2 * np.pi * gA * fit.AXfc
    calibration.kappaAY = 2 * np.pi * gA * fit.AYfc

    calibration.AXoffset = np.mean(X)
    calibration.AYoffset = np.mean(Y)
    calibration.AZoffset = np.mean(Sum)

    fResHigh = 1000
    indRes = np.where((fit.f > fLow) & (fit.f < fResHigh))[0]
    calibration.AXmeanRes = np.mean(fit.AXSpec[indRes] / predictedAX[indRes] - 1)
    calibration.AYmeanRes = np.mean(fit.AYSpec[indRes] / predictedAY[indRes] - 1)

    calibration.AXstdRes = mstd(fit.AXSpec[indRes] / predictedAX[indRes] - 1)
    calibration.AYstdRes = mstd(fit.AYSpec[indRes] / predictedAY[indRes] - 1)

    return calibration, fit
