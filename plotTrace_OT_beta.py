import numpy as np
import matplotlib.pyplot as plt

from matlab_helpers import smooth

# Colours used from the cmu.colors MATLAB package
_CMU_COLORS = {
    'ucla blue': (0.3255, 0.4078, 0.5843),
    'salmon': (0.9804, 0.5020, 0.4471),
    'coral': (1.0000, 0.4980, 0.3137),
}


def _clr(name):
    return np.array(_CMU_COLORS[name])


def _to_float(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return np.nan


def plotTrace_OT_beta(data, ext, stgdata, titletxt):
    """Plot DNA extension (top) and forces / trap positions (bottom) vs time.
    stgdata = (stgTime, stgPos) as returned by StageMoveTime."""
    traceclr = np.vstack([_clr('ucla blue'), _clr('salmon')])
    strDate = data.file[0:6]
    traceNum = data.file[-7:-4]
    stgTime = stgdata[0]
    stgPos = stgdata[1]

    print('no. activity(dip): ' + str(int(np.floor((len(stgTime) - 1) / 4))))

    # %%
    dnaext = ext

    fig = plt.figure(figsize=(16, 9))
    try:
        fig.canvas.manager.set_window_title(f'{strDate}_{traceNum}')
    except Exception:
        pass

    slength = plt.subplot2grid((3, 1), (0, 0), rowspan=2, fig=fig)
    for stgti in range(len(stgTime)):
        pos = _to_float(stgPos[stgti]) if stgti < len(stgPos) else np.nan
        if not np.isnan(pos):  # if stgTime has ended with a stage position
            if pos == 3:
                slength.axvline(stgTime[stgti], color='b', linestyle='--', linewidth=2)
            elif pos == 4:
                slength.axvline(stgTime[stgti], color='r', linestyle='--', linewidth=2)
        else:
            slength.axvline(stgTime[stgti], color='k', linestyle='-', linewidth=1)

    slength.plot(data.time, dnaext, color=traceclr[1, :])
    slength.set_ylabel('DNA Ext (bp)')
    slength.tick_params(labelbottom=False)
    slength.set_title(f'{titletxt} {strDate}_{traceNum}')
    slength.grid(True)

    forces = plt.subplot2grid((3, 1), (2, 0), fig=fig, sharex=slength)
    right = forces.twinx()
    lnwdt1 = 1
    avwin = 666
    right.plot(data.time, data.trappos1, '--', linewidth=lnwdt1, color=_clr('coral'), label='Trappos1')
    right.plot(data.time, data.trappos2, '-', linewidth=lnwdt1, color=_clr('coral'), label='Trappos2')
    right.set_ylabel('AOM(MHz)')

    lnwdt2 = 0.5
    forces.plot(data.time, smooth(data.force_AX, avwin, 'moving'), '-', color=(1, 0, 0), linewidth=lnwdt2, label='Trap AX F')
    forces.plot(data.time, smooth(data.force_BX, avwin, 'moving'), '-', color=(0, 0, 1), linewidth=lnwdt2, label='Trap BX F')
    forces.plot(data.time, data.forceXY, '-', color=(0, 0, 0, 0.1), linewidth=2, label='_nolegend_')
    forces.plot(data.time, smooth(data.forceXY, avwin, 'moving'), '-', color=(0, 0, 0), linewidth=2, label='Avg F')
    forces.plot(data.time, smooth(data.force_AY, avwin, 'moving'), '-.', color=(0, 0.6, 0, 0.5), label='Trap AY F')
    forces.plot(data.time, smooth(data.force_BY, avwin, 'moving'), '-.', color=(0, 0, 0.6, 0.5), label='Trap BY F')

    for stgti in range(len(stgTime)):
        forces.axvline(stgTime[stgti], color='k', linestyle='--')
    forces.set_ylabel('Force (pN)')
    forces.set_xlabel('Time (s)')

    allf = np.concatenate([np.ravel(data.force_AX), np.ravel(data.force_AY),
                           np.ravel(data.force_BX), np.ravel(data.force_BY)])
    miny = np.nanmin(allf)
    maxy = np.nanmax(allf)
    forces.set_ylim([miny - 2, maxy + 2])

    h1, l1 = forces.get_legend_handles_labels()
    h2, l2 = right.get_legend_handles_labels()
    forces.legend(h1 + h2, l1 + l2, fontsize=15)

    return fig
