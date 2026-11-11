"""
CalibrateCont

Looks at today's data folder, and if the newest .dat file has not been
handled yet: runs a calibration (datatype 1), plots a force-extension curve
(datatype 2) or a fluorescence map (datatype 3).

Uses the globals datadirectories, mostrecentcal and calparams
(shared_globals). The MATLAB base-workspace flags (assignin/evalin) are kept
in shared_globals.base_workspace.

Needs Calibrate, plot_force_extension and fplot_fluor_map, which were not
part of the converted folder.
"""
import glob
import os
from datetime import date
import numpy as np
import matplotlib.pyplot as plt

import shared_globals as G


def _f64(fid):
    return float(np.fromfile(fid, dtype='>f8', count=1)[0])


def CalibrateCont():
    D = date.today()

    Date = f'{D.year - 2000:02d}{D.month:02d}{D.day:02d}'

    startpath = os.path.join(G.datadirectories, Date) + os.sep

    F = sorted(glob.glob(startpath + '*.dat'))  # MATLAB dir() order is alphabetical

    if len(F) == 0:
        return

    rootfile = os.path.basename(F[-1])[0:10] + '.dat'
    calfilename = 'cal' + rootfile[0:10]  # this is to be sure to not get confused by the fl files

    ws = G.base_workspace
    if (calfilename not in ws) and (('data_' + rootfile[0:10]) not in ws):
        # try to open the file, will fail if still writing to file
        try:
            with open(startpath + rootfile, 'rb') as fid:
                nheader = _f64(fid)
                version = _f64(fid)       # file version
                rawperiod = _f64(fid)     # raw sampling period of data before averaging for saving (s)
                sampperiod = _f64(fid)    # sampling period of data (s)
                nchannels = _f64(fid)     # number of data channels
                chanpattern = _f64(fid)   # channel pattern identifier (e.g., 1 trap or two)
                datatype = _f64(fid)      # data type identifier (1 = calibration, 0 = time series, 2 = raster scan)
        except Exception:
            return

        fnumber = float(rootfile[7:10])

        if datatype == 1:  # a cal file
            from Calibrate import Calibrate  # not part of the converted folder
            Calibrate(Date, rootfile, G.calparams, 1, 1)
            G.mostrecentcal = float(rootfile[7:10])

        elif datatype == 2:  # a raster scan, use offset, don't run on offset
            if fnumber - G.mostrecentcal > 1:
                try:
                    from plot_force_extension import plot_force_extension  # not part of the converted folder
                    plot_force_extension(Date, G.mostrecentcal, G.mostrecentcal + 1, fnumber)
                    ws['data_' + rootfile[0:10]] = 1
                    plt.gcf().savefig(startpath + rootfile[0:10] + '_fext.png')
                except Exception:
                    ws['data_' + rootfile[0:10]] = 1
                    return

        elif datatype == 3:  # a fluorescence map
            try:
                from fplot_fluor_map import fplot_fluor_map  # not part of the converted folder
                fplot_fluor_map(Date, rootfile)
                plt.gcf().savefig(startpath + rootfile[0:10] + '_image.png')
            except Exception:
                ws['data_' + rootfile[0:10]] = 1
                return

        else:
            ws['data_' + rootfile[0:10]] = 1
