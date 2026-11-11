"""
ReadJointFile

Direct usage example: data = ReadJointFile(rawDataPath + '170328\\', '170328_073.dat')

Read both fleezers and old trap files. Modified from ReadMattFile_v8.
Added ot field to data. Should not affect any fleezer data.
150429 kw

function to read Matt's data files with headers.
090218 mjc

added ability to read trap position files
101025 mjc

All files are big-endian ('ieee-be'), as in the MATLAB version.
`startpath` is joined to `rootfile` by plain string concatenation, exactly
like MATLAB, so startpath should end with a path separator.
"""
import numpy as np
from matlab_helpers import Struct


def _f64(fid):
    return float(np.fromfile(fid, dtype='>f8', count=1)[0])


def ReadJointFile(startpath, rootfile):
    data = Struct()

    # open file
    fid = open(startpath + rootfile, 'rb')

    # read number of header entries following
    _f64(fid)

    # read version number of file format
    data.version = _f64(fid)

    data.rawperiod = _f64(fid)    # raw sampling period of data before averaging for saving (s)
    data.sampperiod = _f64(fid)   # sampling period of data (s)
    data.nchannels = _f64(fid)    # number of data channels
    data.chanpattern = _f64(fid)  # channel pattern identifier
    data.datatype = _f64(fid)     # data type identifier

    # now trap positions, three AOD x,y frequency pairs written, even if only
    # one trap is being used and no interlacing.  These positions are also the
    # starting positions for raster scans.
    data.t1x = _f64(fid)
    data.t1y = _f64(fid)
    data.t2x = _f64(fid)
    data.t2y = _f64(fid)
    data.t3x = _f64(fid)
    data.t3y = _f64(fid)

    # 1 or 3 (odd) if fluorescence data was saved (in another file)
    # 2 or 3 if position data saved
    data.extrasaved = _f64(fid)

    if data.extrasaved == 1:
        data.flsaved = 1
        data.trpossaved = 0
    elif data.extrasaved == 2:
        data.flsaved = 0
        data.trpossaved = 1
    elif data.extrasaved == 3:
        data.flsaved = 1
        data.trpossaved = 1
    else:
        data.flsaved = 0
        data.trpossaved = 0

    if data.version <= 8:
        # next four values are not used yet
        for i in range(4):
            _f64(fid)
        # (as in the original, nothing more is read for version <= 8)
        fid.close()
        return data

    # data.version >= 9
    if data.datatype == 3:  # trap raster scan
        data.scan2dor1d = _f64(fid)
        data.oldtrap = 0
        # next three values are not used yet
        for i in range(3):
            _f64(fid)
    elif data.datatype == 0:  # time series
        data.oldtrap = _f64(fid)
        data.stageX = _f64(fid)
        data.stageY = _f64(fid)
        data.stageZ = _f64(fid)
    else:
        data.oldtrap = _f64(fid)  # 1 if old trap, 0 if fleezers.
        # next three values are not used yet
        for i in range(3):
            _f64(fid)

    # read raster scan parameters if a raster scan
    if data.datatype == 2:  # AOM raster scan
        data.scandir = _f64(fid)    # 1 = x scan, 0 = y scan
        data.scanrange = _f64(fid)  # scan range (MHz), sign gives direction
        data.nsteps = _f64(fid)     # number of frequency steps per scan
        data.ncycles = _f64(fid)    # number of cycles per step
        data.nscanfb = _f64(fid)    # number of scans, forward and back (always come in pairs)
    elif data.datatype == 3:  # MCL raster scan  (AT: movable mirror scans, on the old trap)
        data.scandir = _f64(fid)    # 1 = x scan, 0 = y scan
        data.scanrange = _f64(fid)  # scan range (V), sign gives direction
        data.nsteps = _f64(fid)     # number of voltage steps per scan
        data.ncycles = _f64(fid)    # number of cycles per step
        data.nscanfb = _f64(fid)    # the number of scan lines, forward and back (always come in pairs)
        data.initialx = _f64(fid)   # initial x position (to know where scan starts from)
        data.initialy = _f64(fid)   # initial y position

    # now read a comment string
    nstrbytes = int(_f64(fid))  # number of string bytes
    data.comment = fid.read(nstrbytes).decode('latin-1')

    # now read data
    if data.oldtrap:
        D = np.fromfile(fid, dtype='>f4').astype(float)
    else:
        D = np.fromfile(fid, dtype='>i2').astype(float)

    fid.close()

    # reshape (MATLAB: D = reshape(D, nchannels, [])')
    nch = int(data.nchannels)
    D = D.reshape(-1, nch)

    # convert to voltages
    D = D / 3276.7

    data.path = startpath
    data.file = rootfile

    cp = int(data.chanpattern)
    col = lambda k: D[:, k - 1].copy()  # MATLAB 1-based column

    if cp == 1:
        data.A_Y = col(1)
        data.A_X = col(2)
        data.A_Sum = col(3)

        data.B_Y = col(4)
        data.B_X = col(5)
        data.B_Sum = col(6)

        data.A_FB_Sum = col(7)
        data.B_FB_Sum = col(8)

        data.C_Y = col(9)
        data.C_X = col(10)
        data.C_Y_Sum = col(11)
        data.C_X_Sum = col(12)

        data.time = np.arange(len(data.A_X)) * data.sampperiod

        # normalize position data
        data.A_Y = data.A_Y / data.A_Sum
        data.A_X = data.A_X / data.A_Sum
        data.B_Y = data.B_Y / data.B_Sum
        data.B_X = data.B_X / data.B_Sum
        data.C_Y = data.C_Y / data.C_Y_Sum
        data.C_X = data.C_X / data.C_X_Sum

    elif cp == 2:  # two traps, FB detector X + sum
        data.A_Y = col(1)
        data.A_X = col(2)
        data.A_Sum = col(3)

        data.B_Y = col(4)
        data.B_X = col(5)
        data.B_Sum = col(6)

        data.A_FB_X = -col(7)
        data.A_FB_Sum = -col(8)

        data.B_FB_X = -col(9)
        data.B_FB_Sum = -col(10)

        data.C_Y = col(11)
        data.C_X = col(12)
        data.C_Y_Sum = col(13)
        data.C_X_Sum = -col(13)

        data.time = np.arange(len(data.A_X)) * data.sampperiod

        # normalize position data
        data.A_Y = data.A_Y / data.A_Sum
        data.A_X = data.A_X / data.A_Sum
        data.B_Y = data.B_Y / data.B_Sum
        data.B_X = data.B_X / data.B_Sum
        data.C_Y = data.C_Y / data.C_Y_Sum
        data.C_X = data.C_X / data.C_X_Sum
        data.A_FB_X = data.A_FB_X / data.A_FB_Sum
        data.B_FB_X = data.B_FB_X / data.B_FB_Sum

    elif cp == 4:  # two traps, FB detector X + Y + sum
        data.A_Y = col(1)
        data.A_X = col(2)
        data.A_Sum = col(3)

        data.B_Y = col(4)
        data.B_X = col(5)
        data.B_Sum = col(6)

        data.A_FB_Y = -col(7)
        data.A_FB_X = -col(8)
        data.A_FB_Sum = col(9)

        data.B_FB_Y = -col(10)
        data.B_FB_X = -col(11)
        data.B_FB_Sum = col(12)

        data.C_Y = col(13)
        data.C_X = col(13)
        data.C_Y_Sum = col(14)
        data.C_X_Sum = -col(14)

        data.time = np.arange(len(data.A_X)) * data.sampperiod

        # normalize position data
        data.A_Y = data.A_Y / data.A_Sum
        data.A_X = data.A_X / data.A_Sum
        data.B_Y = data.B_Y / data.B_Sum
        data.B_X = data.B_X / data.B_Sum
        data.C_Y = data.C_Y / data.C_Y_Sum
        data.C_X = data.C_X / data.C_X_Sum
        data.A_FB_Y = data.A_FB_Y / data.A_FB_Sum
        data.A_FB_X = data.A_FB_X / data.A_FB_Sum
        data.B_FB_Y = data.B_FB_Y / data.B_FB_Sum
        data.B_FB_X = data.B_FB_X / data.B_FB_Sum

    elif cp == 3:  # one trap, FB detector position and sum
        data.A_Y = col(1)
        data.A_X = col(2)
        data.A_Sum = col(4)

        # just a copy of A
        data.B_Y = col(1)
        data.B_X = col(2)
        data.B_Sum = col(4)

        data.A_FB_X = -col(8)
        data.A_FB_Sum = col(3)

        data.C_Y = col(5)
        data.C_X = col(6)
        data.C_Y_Sum = col(7)
        data.C_X_Sum = -col(7)

        data.time = np.arange(len(data.A_X)) * data.sampperiod

        # normalize position data
        data.A_Y = data.A_Y / data.A_Sum
        data.A_X = data.A_X / data.A_Sum
        data.B_Y = data.B_Y / data.B_Sum
        data.B_X = data.B_X / data.B_Sum
        data.A_FB_X = data.A_FB_X / data.A_FB_Sum

    elif cp == 0:  # one trap, FB detector sum
        data.A_Y = col(1)
        data.A_X = col(2)
        data.A_Sum = col(4)

        data.A_FB_Sum = col(3)

        data.B_Y = col(5)
        data.B_X = col(6)
        data.B_Y_Sum = col(7)
        data.B_X_Sum = col(8)
        data.B_Sum = col(7)

        data.time = np.arange(len(data.A_X)) * data.sampperiod

        # normalize position data
        data.A_Y = data.A_Y / data.A_Sum
        data.A_X = data.A_X / data.A_Sum
        data.B_Y = data.B_Y / data.B_Sum
        data.B_X = data.B_X / data.B_Sum

    elif cp == 5:  # old trap data
        data.A_Y = col(1)
        data.A_X = col(2)
        data.A_Sum = col(5)

        data.B_Y = col(3)
        data.B_X = col(4)
        data.B_Sum = col(6)

        data.time = np.arange(len(data.A_X)) * data.sampperiod

        # normalize position data
        data.A_Y = data.A_Y / data.A_Sum
        data.A_X = data.A_X / data.A_Sum
        data.B_Y = data.B_Y / data.B_Sum
        data.B_X = data.B_X / data.B_Sum

    # now read fluorescence data file if it exists
    if data.flsaved == 1:
        fid = open(startpath + rootfile[:-4] + '_fl.dat', 'rb')

        # read number of header entries following
        _f64(fid)

        data.flintperiod = _f64(fid)    # APD integration period of data (s)
        data.napdsampstep = _f64(fid)   # number of APD samples per scan step (only meaningful for a scan)

        # next 6 entries are reserved for now
        for i in range(6):
            _f64(fid)

        # now read data to the end of the file, the two APD's are interlaced
        F = np.fromfile(fid, dtype='>u4').astype(float)
        fid.close()

        # reshape to separate APD 1 and 2
        F = F.reshape(-1, 2)

        data.apd1 = F[:, 0].copy()
        data.apd2 = F[:, 1].copy()

        data.apdtime = np.arange(1, len(data.apd1) + 1) * data.flintperiod

    # now read trap position data file if it exists
    if data.trpossaved == 1:
        fid = open(startpath + rootfile[:-4] + '_pos.dat', 'rb')

        # read number of header entries following
        _f64(fid)

        # next 8 entries are reserved for now
        for i in range(8):
            _f64(fid)

        # now read data to the end of the file, there are two interlaced trap
        # positions (assumed for now)
        if data.oldtrap:
            F = np.fromfile(fid, dtype='>f4').astype(float)
        else:
            F = np.fromfile(fid, dtype='>u8').astype(float)  # might work

        # reshape to separate traps 1 and 2,
        if data.oldtrap:
            # 191010 SY added. not for using for gen pop. if problem let me know.
            if len(F) == len(data.A_X) * 4:
                F = F.reshape(-1, 4)
                data.trappos1 = F[:, 0].copy()
                data.trappos2 = F[:, 1].copy()
                data.offsetfitY = F[:, 2].copy()
                data.offsetfitX = F[:, 3].copy()
            elif len(F) == len(data.A_X) * 5:
                F = F.reshape(-1, 5)
                data.trappos1 = F[:, 0].copy()
                data.trappos2 = F[:, 1].copy()
                data.offsetfitX = F[:, 2].copy()
                data.offsetfitY = F[:, 3].copy()
                data.imgsumsig = F[:, 4].copy()
            else:
                F = F.reshape(-1, 2)
                data.trappos1 = F[:, 0].copy()
                data.trappos2 = F[:, 1].copy()
        else:
            F = F.reshape(-1, 2)
            data.trappos1 = F[:, 0] * 49.152e6 * 6 / (1e6 * 2**48)  # convert to frequency
            data.trappos2 = F[:, 1] * 49.152e6 * 6 / (1e6 * 2**48)
        fid.close()

    return data
