import copy
import numpy as np
from matlab_helpers import Struct, isfield, lfilter_ma, interp1


def force_offset_subtraction_only_Data_SY(offsetData, trapData, calData):
    """Subtract the QPD force offset (measured in offsetData) from trapData in
    both X and Y, then convert to nm, DNA extension and force."""
    trapData = copy.copy(trapData)  # MATLAB passes structs by value
    trapData.calData = calData

    # %% testing out offset subtraction in both X Y
    TrapSepOff = np.asarray(offsetData.TrapSep, dtype=float).ravel()
    sid = np.argsort(TrapSepOff, kind='stable')
    s_extOffset = TrapSepOff[sid]
    s_Offset = Struct()
    s_Offset.A_X = np.asarray(offsetData.A_X).ravel()[sid]
    s_Offset.B_X = np.asarray(offsetData.B_X).ravel()[sid]
    s_Offset.A_Y = np.asarray(offsetData.A_Y).ravel()[sid]
    s_Offset.B_Y = np.asarray(offsetData.B_Y).ravel()[sid]

    # Smooth offset
    Nfilt = int(round(1 / (trapData.sampperiod * 10)))  # Divide by 10 to avoid oversmoothing

    ss_extOffset_X = lfilter_ma(Nfilt, s_extOffset)
    ss_Offset = Struct()
    ss_Offset.A_X = lfilter_ma(Nfilt, s_Offset.A_X)
    ss_Offset.B_X = lfilter_ma(Nfilt, s_Offset.B_X)

    ss_extOffset_Y = lfilter_ma(Nfilt, s_extOffset)
    ss_Offset.A_Y = lfilter_ma(Nfilt, s_Offset.A_Y)
    ss_Offset.B_Y = lfilter_ma(Nfilt, s_Offset.B_Y)

    # Filter offset (MATLAB x(Nfilt:end) -> python x[Nfilt-1:])
    k = Nfilt - 1
    sss_extOffset_X = ss_extOffset_X[k:]
    sss_extOffset_Y = ss_extOffset_Y[k:]
    sss_Offset = Struct()
    sss_Offset.A_X = ss_Offset.A_X[k:]
    sss_Offset.A_Y = ss_Offset.A_Y[k:]
    sss_Offset.B_X = ss_Offset.B_X[k:]
    sss_Offset.B_Y = ss_Offset.B_Y[k:]

    ssss_extOffset_X, index = np.unique(sss_extOffset_X, return_index=True)
    sss_Offset.A_X = sss_Offset.A_X[index]
    sss_Offset.B_X = sss_Offset.B_X[index]

    ssss_extOffset_Y, index = np.unique(sss_extOffset_Y, return_index=True)
    sss_Offset.A_Y = sss_Offset.A_Y[index]
    sss_Offset.B_Y = sss_Offset.B_Y[index]

    # Interpolate offset values corresponding to data values
    ssi_Offset = Struct()
    ssi_Offset.A_X = interp1(ssss_extOffset_X, sss_Offset.A_X, trapData.TrapSep)
    ssi_Offset.B_X = interp1(ssss_extOffset_X, sss_Offset.B_X, trapData.TrapSep)

    ssi_Offset.A_Y = interp1(ssss_extOffset_Y, sss_Offset.A_Y, trapData.TrapSep)
    ssi_Offset.B_Y = interp1(ssss_extOffset_Y, sss_Offset.B_Y, trapData.TrapSep)

    # Replace the raw values with force offset subtracted values
    trapData.A_X = trapData.A_X - ssi_Offset.A_X
    trapData.B_X = trapData.B_X - ssi_Offset.B_X
    trapData.A_Y = trapData.A_Y - ssi_Offset.A_Y
    trapData.B_Y = trapData.B_Y - ssi_Offset.B_Y

    # %% CONVERT QPD DATA FROM VOLTS TO NANOMETERS
    trapData.A_Xnm = (trapData.A_X) * calData.alphaAX
    trapData.B_Xnm = (-trapData.B_X) * calData.alphaBX
    trapData.A_Ynm = (-trapData.A_Y) * calData.alphaAY
    trapData.B_Ynm = (trapData.B_Y) * calData.alphaBY

    # %% DNA EXTENSION IN NANOMETERS
    trapData.DNAext_X = trapData.TrapSep - trapData.A_Xnm - trapData.B_Xnm - calData.beadA / 2 - calData.beadB / 2
    trapData.DNAext_Y = trapData.TrapSep - trapData.A_Ynm - trapData.B_Ynm - calData.beadA / 2 - calData.beadB / 2

    if isfield(trapData, 'TrapSepXY'):
        trapData.BtoB_X = trapData.TrapSep_Xnm - (trapData.A_Xnm + trapData.B_Xnm)
        trapData.BtoB_Y = trapData.TrapSep_Ynm - (trapData.A_Ynm + trapData.B_Ynm)
        trapData.DNAext = np.sqrt(trapData.BtoB_X**2 + trapData.BtoB_Y**2) - calData.beadA / 2 - calData.beadB / 2

    # %% FORCE (pN)
    trapData.force_AX = trapData.A_Xnm * calData.kappaAX
    trapData.force_BX = trapData.B_Xnm * calData.kappaBX
    trapData.force_AY = trapData.A_Ynm * calData.kappaAY
    trapData.force_BY = trapData.B_Ynm * calData.kappaBY
    trapData.forceX = (trapData.force_AX + trapData.force_BX) / 2
    trapData.forceY = (trapData.force_AY + trapData.force_BY) / 2
    trapData.forceXY = np.sqrt(trapData.forceX**2 + trapData.forceY**2)

    return trapData
