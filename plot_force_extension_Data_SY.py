import copy
import numpy as np
import matplotlib.pyplot as plt

import shared_globals as G
from matlab_helpers import Struct, isfield, movmean, logical_index, findpeaks
from ConstructConstants import ConstructConstants
from XWLCContour import XWLCContour
from getFakeoffset import getFakeoffset


def plot_force_extension_Data_SY(Data, construct=None):
    """
    Plot force vs. extension for a raster-scan (pulling) trace together with
    the WLC model(s) of the construct.

    Uses the global WLC_param. In MATLAB `construct` was also declared global;
    here the argument is used, falling back to shared_globals.construct.
    Returns outData (a copy of Data with ConstructConst, meanFpull and
    meanFrelax added).
    """
    if construct is None:
        construct = G.construct
    WLC_param = G.WLC_param

    outData = copy.copy(Data)
    Data = copy.copy(Data)  # MATLAB passes structs by value
    calData = Data.calData

    param = Struct()
    param.Date = Data.file[0:6]
    param.calfilenumber = Data.rootCalFile
    param.fextnumber = Data.rootTraceFile
    if Data.scandir == 0:
        param.direction = 'Y'  # Y or X
    elif Data.scandir == 1:
        param.direction = 'X'  # Y or X

    # %% load construct parameters
    ConstructConst = ConstructConstants(construct)
    print('1')
    outData.ConstructConst = ConstructConst

    # %% plot the WLC model
    max_force = np.nanmax(np.concatenate([np.ravel(Data.force_AX), np.ravel(Data.force_AY),
                                       np.ravel(Data.force_BX), np.ravel(Data.force_BY)]))
    F = np.arange(0, max_force + 1e-12, 0.05)
    bp = 1 * WLC_param.hds * XWLCContour(F, WLC_param.Pds, WLC_param.Sds, WLC_param.kT)  # nm/bp
    nt = 1 * WLC_param.hss * XWLCContour(F, WLC_param.Pss, WLC_param.Sss, WLC_param.kT)  # nm/nt

    fig = plt.figure(figsize=(10.5, 6.59))
    try:
        fig.canvas.manager.set_window_title(param.fextnumber[:-4])
    except Exception:
        pass
    ax = fig.gca()
    if ConstructConst.hairpin == 1:
        # closed model
        WLC_ext = ConstructConst.ssDNA1 * nt + ConstructConst.dsDNA1 * bp
        # open model
        WLC_ext_op = ConstructConst.ssDNA2 * nt + ConstructConst.dsDNA2 * bp

        ax.plot(WLC_ext, F, '--r', linewidth=2, label='WLC closed')
        ax.plot(WLC_ext_op, F, '--k', linewidth=2, label='WLC open')
    else:
        WLC_ext = ConstructConst.ssDNA1 * nt + ConstructConst.dsDNA1 * bp
        ax.plot(WLC_ext, F, '--r', linewidth=2, label='model')

    # %% plot the data

    # Filter TrapSep
    TrapSepFilt = movmean(Data.TrapSep, 30)

    if isfield(Data, 'TrapSepXY'):
        TrapSepFilt = movmean(Data.TrapSepXY, 30)

    # remove jumps in trap position
    TrapSepFilt_diff = np.diff(TrapSepFilt)
    normalized_diff = np.diff(TrapSepFilt) / np.mean(TrapSepFilt_diff[TrapSepFilt_diff > 0])
    keep = (normalized_diff > 0.8) | (normalized_diff < -0.8)
    TrapSepFiltClean = logical_index(TrapSepFilt, keep)
    for name in ['force_AX', 'force_BX', 'force_AY', 'force_BY', 'DNAext_X', 'DNAext_Y',
                 'BtoB_X', 'BtoB_Y', 'DNAext', 'forceXY', 'forceX', 'forceY']:
        if isfield(Data, name):  # BtoB_*/DNAext only exist for old-trap XY data
            Data[name] = logical_index(Data[name], keep)

    # find the top peaks
    numMAX, indMAX = findpeaks(TrapSepFiltClean)
    # find the bottom peaks
    TrapSepInverse = np.concatenate([[-np.max(numMAX)], -TrapSepFiltClean, [-np.max(numMAX)]])
    numMIN, indMIN = findpeaks(TrapSepInverse)
    indMIN = indMIN - 1
    # count number of raster scans
    rasterNUM = len(indMAX)

    # %% get the force array
    if param.direction == 'X':
        FA = Data.force_AX
        FB = Data.force_BX
        FM = (Data.force_AX + Data.force_BX) / 2
        DNAext = Data.DNAext_X
    elif param.direction == 'Y':
        FA = Data.force_AY
        FB = Data.force_BY
        FM = (Data.force_AY + Data.force_BY) / 2
        DNAext = Data.DNAext_Y

    if isfield(Data, 'TrapSepXY'):
        DNAext = Data.DNAext
        FM = Data.forceXY
        FA = np.sqrt(Data.force_AX**2 + Data.force_AY**2)
        FB = np.sqrt(Data.force_BX**2 + Data.force_BY**2)

    # %% get fake offset
    # get usual fakeoffset to get initial fminsearch value
    fr = ConstructConst.fitrange
    DNAext_offset = np.mean(DNAext[(FM > fr[0]) & (FM < fr[1])])
    WLC_offset = np.mean(WLC_ext[(F > fr[0]) & (F < fr[1])])
    fake_offset_0 = DNAext_offset - WLC_offset
    fake_offset = getFakeoffset(DNAext, FM, fr, fake_offset_0, construct)
    DNAext = DNAext - fake_offset

    # %%
    orange = [0.9100, 0.4100, 0.1700]
    blue = [0, 0.4470, 0.7410]
    outData.meanFpull = [None] * rasterNUM
    outData.meanFrelax = [None] * rasterNUM
    for ii in range(rasterNUM):
        indxS = indMIN[ii]
        indxE = indMAX[ii]
        indxEE = indMIN[ii + 1]
        pull = slice(indxS, indxE + 1)     # MATLAB indxS:indxE (inclusive)
        relax = slice(indxE, indxEE + 1)   # MATLAB indxE:indxEE
        ax.plot(DNAext[pull], FM[pull], color=orange, linewidth=1.5, label='mean F pull')
        ax.plot(DNAext[relax], FM[relax], color=blue, linewidth=1.5, label='mean F relax')

        ax.plot(DNAext[pull], FA[pull], '-.', color=orange, linewidth=0.5, label='FA pull')
        ax.plot(DNAext[relax], FA[relax], '-.', color=blue, linewidth=0.5, label='FA relax')
        ax.plot(DNAext[pull], FB[pull], '-.', color=orange, linewidth=0.5, label='FB pull')
        ax.plot(DNAext[relax], FB[relax], '-.', color=blue, linewidth=0.5, label='FB relax')

        outData.meanFpull[ii] = np.column_stack([np.ravel(DNAext[pull]), np.ravel(FM[pull])])
        outData.meanFrelax[ii] = np.column_stack([np.ravel(DNAext[relax]), np.ravel(FM[relax])])

    ax.set_xlim([np.nanmin(DNAext) - 100, np.nanmax(DNAext) + 100])  # MATLAB min/max ignore NaN
    ax.set_ylim([min(0, np.nanmin(FM) - 5), np.nanmax(FM) + 5])
    ax.set_xlabel('DNA extension (nm)')
    ax.set_ylabel('Force (pN)')
    if Data.oldtrap == 1:
        instrument = 'Old-Trap'
    elif Data.oldtrap == 0:
        instrument = 'Fleezers'
    if rasterNUM == 1:
        if ConstructConst.hairpin == 1:
            labels = ['Hairpin closed model', 'Hairpin open model', 'DNA ext. pull, mean force',
                      'DNA ext. relax, mean force', 'DNA ext. pull, A bead force', 'DNA ext. relax, A bead force',
                      'DNA ext. pull, B bead force', 'DNA ext. relax, B bead force']
        else:
            labels = ['WLC model', 'DNA ext. pull, mean force', 'DNA ext. relax, mean force',
                      'DNA ext. pull, A bead force', 'DNA ext. relax, A bead force',
                      'DNA ext. pull, B bead force', 'DNA ext. relax, B bead force']
        ax.legend(ax.lines[:len(labels)], labels, loc='upper left')

    ax.set_title(f'Force vs. Extension, {fake_offset:g} nm offset ({param.fextnumber})\n'
                 f'{construct} construct, {instrument}')

    return outData
