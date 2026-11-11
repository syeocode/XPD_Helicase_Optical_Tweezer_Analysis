"""
This code processes trap data.
02-04-2018 original
DO NOT MODIFY!!
Updated 190520 SY
take in data not file numbers

Input variables:
    calData: data structure output from Calibrate_Data_SY
    traceData: data structure from ReadJointFile
    save: 1 to plot (and save) an overview figure
    beta: 2x2 V-to-pixel matrix (None or [] for no mixing)
Output structure array called traceData contains:
    traceData.TrapSep: trap separation in nm (used conversion factors - see below)
    traceData.A_X / B_X / A_Y / B_Y: QPD data in volts
    traceData.A_Xnm / B_Xnm / A_Ynm / B_Ynm: QPD data in nanometers
    traceData.DNAext_X / DNAext_Y: DNA extension in nanometers in x / y
    traceData.force_AX / force_BX / force_AY / force_BY: force in pN
    traceData.time: time in seconds
Note: To get trap separation on the old-trap, we use Pythagoras's Theorem.
The AOM on the fleezers is setup in X only.
    traceData.trappos1 - moving trap position in X
    traceData.trappos2 - moving trap position in Y
    traceData.t1x / t1y / t2x / t2y: start positions
Conversion factors from voltage or frequency to nm, as of 02-04-2018:
    flconv = 123: nm/MHz from grid measurements for fleezers
    otconvX = 891.5: nm/V from grid measurements for old-trap in X
    otconvY = 505.7: nm/V from grid measurements for old-trap in Y
"""
import copy
import os
import numpy as np
import matplotlib.pyplot as plt

import shared_globals as G
from matlab_helpers import isfield, isempty
from getConvFactor import getConvFactor


def trap_data_Data_sy_2023(calData, traceData, save=0, beta=None):
    traceData = copy.copy(traceData)  # MATLAB passes structs by value
    strDate = traceData.file[0:6]
    traceData.rootTraceFile = traceData.file
    traceData.rootCalFile = calData.file

    ConvFactor = getConvFactor(strDate)

    flconv = ConvFactor.flconv    # nm/MHz
    otconvX = ConvFactor.otconvX  # nm/V 230810, 180629: 883.59;
    otconvY = ConvFactor.otconvY  # nm/V 230810, 180629: 508.93;

    if isfield(traceData, 'scandir'):
        if traceData.scandir == 1:
            traceData.direction = 'X'
        elif traceData.scandir == 0:
            traceData.direction = 'Y'

    # %% trap position
    traceData.trappos1 = np.asarray(traceData.trappos1, dtype=float).ravel()
    traceData.trappos2 = np.asarray(traceData.trappos2, dtype=float).ravel()

    # %% V or MHz CONVERSION TO nm
    if traceData.oldtrap == 0:  # FL
        traceData.TrapSep = (traceData.trappos2 - traceData.trappos1) * flconv
        traceData.flconv = flconv
    elif traceData.oldtrap == 1:  # OT
        # new code 230829 to convert from tilted sm
        traceData.otconvX = otconvX
        traceData.otconvY = otconvY  # keep record of conv factor
        TrapSep_Vx = (traceData.trappos1 - traceData.t2x).reshape(1, -1)
        TrapSep_Vy = (traceData.trappos2 - traceData.t2y).reshape(1, -1)
        TrapSep_V = np.vstack([TrapSep_Vx, TrapSep_Vy])  # if sm tilted, sm V to pixel frame will include both x and y comp

        if isempty(beta):
            # if no beta is inputed use beta with no mixing
            beta = np.array([[1, 0],
                             [0, 1]], dtype=float)
        beta = np.asarray(beta, dtype=float)

        TrapSep_P = beta @ TrapSep_V  # convert V to P using beta

        pixperV_X = np.linalg.norm(beta @ np.array([1, 0]))
        pixperV_Y = np.linalg.norm(beta @ np.array([0, 1]))  # get 1V step in x and y to get pixel per volts
        nmperpix_X = otconvX / pixperV_X
        nmperpix_Y = otconvY / pixperV_Y  # nm/px. oTconv = nm/px( from cal grid ) * px/V ( from bead grid ).
        traceData.TrapSep_Xnm = np.abs(TrapSep_P[0, :] * nmperpix_X)  # Trapsep V to nm
        traceData.TrapSep_Ynm = np.abs(TrapSep_P[1, :] * nmperpix_Y)  # Trapsep V to nm
        traceData.TrapSepXY = np.sqrt(traceData.TrapSep_Xnm**2 + traceData.TrapSep_Ynm**2)
        traceData.TrapSep = np.sqrt(((traceData.trappos1 - traceData.t2x) * otconvX)**2 +
                                    ((traceData.trappos2 - traceData.t2y) * otconvY)**2)

    # %% QPD OFFSET SUBTRACTION
    traceData.A_X = np.asarray(traceData.A_X) - calData.AXoffset  # both are in volts
    traceData.B_X = np.asarray(traceData.B_X) - calData.BXoffset
    traceData.A_Y = np.asarray(traceData.A_Y) - calData.AYoffset
    traceData.B_Y = np.asarray(traceData.B_Y) - calData.BYoffset

    # %% QPD DATA IN NANOMETERS
    traceData.A_Xnm = (traceData.A_X) * calData.alphaAX
    traceData.B_Xnm = (-traceData.B_X) * calData.alphaBX
    traceData.A_Ynm = (-traceData.A_Y) * calData.alphaAY
    traceData.B_Ynm = (traceData.B_Y) * calData.alphaBY

    # %% DNA EXTENSION IN NANOMETERS
    traceData.DNAext_X = traceData.TrapSep - traceData.A_Xnm - traceData.B_Xnm - calData.beadA / 2 - calData.beadB / 2
    traceData.DNAext_Y = traceData.TrapSep - traceData.A_Ynm - traceData.B_Ynm - calData.beadA / 2 - calData.beadB / 2

    if isfield(traceData, 'TrapSepXY'):
        traceData.BtoB_X = traceData.TrapSep_Xnm - (traceData.A_Xnm + traceData.B_Xnm)
        traceData.BtoB_Y = traceData.TrapSep_Ynm - (traceData.A_Ynm + traceData.B_Ynm)
        traceData.DNAext = np.sqrt(traceData.BtoB_X**2 + traceData.BtoB_Y**2) - calData.beadA / 2 - calData.beadB / 2

    # %% FORCE (pN)
    traceData.force_AX = traceData.A_Xnm * calData.kappaAX
    traceData.force_BX = traceData.B_Xnm * calData.kappaBX
    traceData.force_AY = traceData.A_Ynm * calData.kappaAY
    traceData.force_BY = traceData.B_Ynm * calData.kappaBY
    traceData.forceX = (traceData.force_AX + traceData.force_BX) / 2
    traceData.forceY = (traceData.force_AY + traceData.force_BY) / 2
    traceData.forceXY = np.sqrt(traceData.forceX**2 + traceData.forceY**2)

    # %% PLOT
    if save == 1:
        # (the MATLAB version referred to undefined Date/datafnumber/type here;
        #  the trace file name is used instead)
        name = traceData.file[:-4]
        fig, axs = plt.subplots(3, 1, figsize=(10, 9))
        axs[0].plot(traceData.time, traceData.TrapSep)
        axs[0].set_ylabel('Trap Separation (nm)')
        axs[0].legend(['Trap Separation'])
        axs[0].set_title(f'Trap Position vs. Time ({name})')
        axs[1].plot(traceData.time, traceData.DNAext_X)
        axs[1].plot(traceData.time, traceData.DNAext_Y)
        axs[1].set_ylabel('DNA extension (nm)')
        axs[1].legend(['DNA extension in x', 'DNA extension in y'])
        axs[2].plot(traceData.time, traceData.force_AX)
        axs[2].plot(traceData.time, traceData.force_BX)
        axs[2].plot(traceData.time, traceData.force_AY)
        axs[2].plot(traceData.time, traceData.force_BY)
        axs[2].set_ylabel('Force (pN)')
        axs[2].legend(['force bead A in x', 'force bead B in x', 'force bead A in y', 'force bead B in y'])
        fig.savefig(os.path.join(G.savedirectory, name + '.png'))

    return traceData
