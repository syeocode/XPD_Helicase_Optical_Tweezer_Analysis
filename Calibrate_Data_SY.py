"""
This is the shared driver code for doing calibrations using either the
old trap or the fleezers. DO NOT MODIFY THIS CODE!
240102 Modified code to take data SY
This function returns a structure variable 'calData' which contains:
    - bead: bead diameters for bead A and B
    - alpha: the V to nm conversion factor for each bead in X and Y
    - kappa: trap stiffness for each bead in X and Y
    - offset: QPD offset for each bead in X and in Y
The input variables are:
    - data: data structure from ReadJointFile
    - calarray: the parameters to get the calibration file, e.g. calparams = [900 880 125000 20000 8000 1200 100]
    - showplots: an option to plot the calibration file, 0 - do not plot, 1 - plot
    - writecal: an option to save the calibration file, 0 - no not save, 1 - save
"""
import numpy as np
import matplotlib.pyplot as plt

import shared_globals as G
from matlab_helpers import Struct
from getnumbers import getnumbers
from CalibrateExact_Stripped_Data_SY import CalibrateExact_Stripped_Data_SY


def Calibrate_Data_SY(data, calarray=None, showplots=1, writecal=1):
    calparams = G.calparams

    if calarray is None:
        prompts = ['Bead A diameter (nm)', 'Bead B diameter (nm)', 'f_{samp} (Hz)',
                   'f_{XY High} (Hz)', 'f_{Sum High} (Hz)', 'f_{Low} (Hz)', 'Averaging window']
        defaults = [str(v) for v in calparams[:7]]
        calarray = getnumbers('Enter parameters:', prompts, defaults)

    calData = Struct()
    calData.file = data.file
    calData.path = data.path

    # Change array of input parameters to structure variable (better for saving
    # user data later)
    calpar = Struct()
    calpar.beadA = calarray[0]
    calpar.beadB = calarray[1]
    calpar.f_samp = calarray[2]
    calpar.f_xyhi = calarray[3]
    calpar.f_sumhi = calarray[4]
    calpar.f_low = calarray[5]
    calpar.avwin = calarray[6]

    # which signal sets to analyze
    calsets = [1,   # Trap A
               1,   # Trap B
               0]   # Detection laser

    param = np.zeros(7)
    param[6] = calpar.avwin    # ave windows
    param[4] = calpar.f_sumhi  # f sum high

    allfit = Struct()

    if calsets[0] == 1:  # get cal for trap 1 (A)
        param[0] = calpar.beadA
        param[3] = calpar.f_xyhi  # fxy high
        param[5] = calpar.f_low   # f low
        fitmodel = 'AliasedHydro'
        data_sub = Struct(sampperiod=data.sampperiod, X=data.A_X, Y=data.A_Y,
                          Sum=data.A_Sum, path=data.path, file=data.file)

        cal, fit = CalibrateExact_Stripped_Data_SY(param, fitmodel, data_sub)

        calData.beadA = param[0]
        calData.alphaAX = cal.alphaAX
        calData.alphaAY = cal.alphaAY
        calData.kappaAX = cal.kappaAX
        calData.kappaAY = cal.kappaAY
        calData.AXoffset = cal.AXoffset
        calData.AYoffset = cal.AYoffset

        allfit.fA = fit.f
        allfit.AXSpecRaw = fit.AXSpecRaw
        allfit.AYSpecRaw = fit.AYSpecRaw
        allfit.predictedAX = fit.predictedAX
        allfit.predictedAY = fit.predictedAY

    if calsets[1] == 1:  # get cal for trap 2 (B)
        param[0] = calpar.beadB
        param[3] = calpar.f_xyhi  # fxy high
        param[5] = calpar.f_low   # f low
        fitmodel = 'AliasedHydro'
        data_sub = Struct(sampperiod=data.sampperiod, X=data.B_X, Y=data.B_Y,
                          Sum=data.B_Sum, path=data.path, file=data.file)

        cal, fit = CalibrateExact_Stripped_Data_SY(param, fitmodel, data_sub)

        # All "cal" fields are named as if they were bead A from the
        # "CalibrateExact_Stripped" function, but they need to be renamed here,
        # since this is actually the calibration for bead B.
        calData.beadB = param[0]
        calData.alphaBX = cal.alphaAX
        calData.alphaBY = cal.alphaAY
        calData.kappaBX = cal.kappaAX
        calData.kappaBY = cal.kappaAY
        calData.BXoffset = cal.AXoffset
        calData.BYoffset = cal.AYoffset

        allfit.fB = fit.f
        allfit.BXSpecRaw = fit.AXSpecRaw
        allfit.BYSpecRaw = fit.AYSpecRaw
        allfit.predictedBX = fit.predictedAX
        allfit.predictedBY = fit.predictedAY

    if calsets[2] == 1:  # get cal for detection laser on trap
        param[0] = 1000   # Bead B
        param[3] = 15e3   # fxy high
        param[5] = 200    # f low
        fitmodel = 'AliasedFiltered'
        data_sub = Struct(sampperiod=data.sampperiod, X=data.C_X, Y=data.C_Y,
                          Sum=data.C_Y_Sum, path=data.path, file=data.file)

        cal, fit = CalibrateExact_Stripped_Data_SY(param, fitmodel, data_sub)

        calData.beadC = param[0]
        calData.alphaCX = cal.alphaAX
        calData.alphaCY = cal.alphaAY
        calData.kappaCX = cal.kappaAX
        calData.kappaCY = cal.kappaAY
        calData.CXoffset = cal.AXoffset
        calData.CYoffset = cal.AYoffset

        allfit.fC = fit.f
        allfit.CXSpecRaw = fit.AXSpecRaw
        allfit.CYSpecRaw = fit.AYSpecRaw
        allfit.predictedCX = fit.predictedAX
        allfit.predictedCY = fit.predictedAY

    if showplots:
        columns = sum(calsets)
        rootfile = calData.file
        fig, axs = plt.subplots(2, columns, squeeze=False, figsize=(5 * columns, 8))
        try:
            fig.canvas.manager.set_window_title('Calibration for ' + rootfile[:-4])
        except Exception:
            pass
        prev = 0
        if calsets[0] == 1:
            ax = axs[0, 0]
            ax.loglog(allfit.fA, allfit.AXSpecRaw)
            ax.loglog(allfit.fA, allfit.predictedAX, 'k')
            ax.set_title(rf'AX: $\kappa$ = {round(calData.kappaAX, 3):g}, $\alpha$ = '
                         rf'{round(calData.alphaAX, 2):g}, Offset = {round(calData.AXoffset, 2):g}')
            ax.set_ylabel(r'Power (V$^2$ s)')
            ax = axs[1, 0]
            ax.loglog(allfit.fA, allfit.AYSpecRaw)
            ax.loglog(allfit.fA, allfit.predictedAY, 'k')
            ax.set_title(rf'AY: $\kappa$ = {round(calData.kappaAY, 2):g}, $\alpha$ = '
                         rf'{round(calData.alphaAY, 2):g}, Offset = {round(calData.AYoffset, 2):g}')
            prev = 1
            ax.set_xlabel('Frequency (Hz)')
            ax.set_ylabel(r'Power (V$^2$ s)')

        if calsets[1] == 1:
            ax = axs[0, prev]
            ax.loglog(allfit.fB, allfit.BXSpecRaw)
            ax.loglog(allfit.fB, allfit.predictedBX, 'k')
            ax.set_title(rf'BX: $\kappa$ = {round(calData.kappaBX, 2):g}, $\alpha$ = '
                         rf'{round(calData.alphaBX, 2):g}, Offset = {round(calData.BXoffset, 2):g}')
            ax = axs[1, prev]
            ax.loglog(allfit.fB, allfit.BYSpecRaw)
            ax.loglog(allfit.fB, allfit.predictedBY, 'k')
            ax.set_title(rf'BY: $\kappa$ = {round(calData.kappaBY, 2):g}, $\alpha$ = '
                         rf'{round(calData.alphaBY, 2):g}, Offset = {round(calData.BYoffset, 2):g}')
            prev = prev + 1
            ax.set_xlabel('Frequency (Hz)')

        if calsets[2] == 1:
            ax = axs[0, prev]
            ax.loglog(allfit.fC, allfit.CXSpecRaw)
            ax.loglog(allfit.fC, allfit.predictedCX, 'k')
            ax.set_title(rf'CX: $\kappa$ = {round(calData.kappaCX, 2):g}, $\alpha$ = '
                         rf'{round(calData.alphaCX, 2):g}, Offset = {round(calData.CXoffset, 2):g}')
            ax = axs[1, prev]
            ax.loglog(allfit.fC, allfit.CYSpecRaw)
            ax.loglog(allfit.fC, allfit.predictedCY, 'k')
            ax.set_title(rf'CY: $\kappa$ = {calData.kappaCY:.3g}. $\alpha$ = '
                         rf'{calData.alphaCY:.4g}, Offset = {calData.CYoffset:.3g}')
        fig.tight_layout()

    return calData
