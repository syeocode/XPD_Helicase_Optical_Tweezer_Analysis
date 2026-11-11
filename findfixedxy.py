"""
findfixedxy(date)

Find the fixed-trap bead position (in pixels and in movable-trap volts) from
bead-grid images, and the linear transformation `beta` (V -> pixel).
Results are saved to <datadirectories><date>\\fixedxyres_<date>.mat.

Needs `beadcentercalc` (not part of the converted folder) importable as
    from beadcentercalc import beadcentercalc
returning an object with fields x0, y0, pcorr, ssep, r2.
"""
import glob
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from scipy.interpolate import griddata
from scipy.io import savemat

import shared_globals as G
from matlab_helpers import mstd, inputdlg

OUTPUT_FIXED_FIELD = 'fixed'  # name of the fixed-position field in `output`


def _im2double(im):
    im = np.asarray(im)
    if im.dtype == np.uint8:
        return im.astype(float) / 255.0
    if im.dtype == np.uint16:
        return im.astype(float) / 65535.0
    if im.dtype == bool:
        return im.astype(float)
    return im.astype(float)


def _imread(path):
    return np.asarray(Image.open(path))


def _put_slice(x, i, im):
    """MATLAB x(:,:,i) = im (i is 0-based here); grows the stack if needed."""
    if x is None or x.size == 0:
        x = np.zeros(im.shape + (0,))
    if i >= x.shape[2]:
        x = np.concatenate([x, np.zeros(im.shape + (i + 1 - x.shape[2],))], axis=2)
    x[:, :, i] = im
    return x


def findfixedxy(date):
    from beadcentercalc import beadcentercalc  # external dependency

    datadirectories = G.datadirectories

    date = str(date)
    strdate = date
    imdir = os.path.join(datadirectories, date) + os.sep
    print(imdir)
    extra_title = ' fixed x/y based on grid positions'
    imagesavename = date + '_fixedxy'

    # %% get background
    imamp = 1  # 1 or 255
    x = None
    files = sorted(glob.glob(os.path.join(imdir, '*background*.bmp')))
    for i, fn in enumerate(files):
        im = imamp * _im2double(_imread(fn))
        x = _put_slice(x, i, im)
    bgp = np.mean(x, axis=2)
    print('background done')

    # %% get moving images
    xarr = []
    yarr = []
    files = sorted(glob.glob(os.path.join(imdir, '*GridMove*.bmp')))
    for fn in files:
        strarr = os.path.basename(fn).split('_')
        xarr.append(strarr[2])
        yarr.append(strarr[3])

    xarr = sorted(set(xarr))
    yarr = sorted(set(yarr))

    # %% get fixed trap bead position
    print('fixed trap calc.. ')
    x = None
    files = sorted(glob.glob(os.path.join(imdir, '*fixed*.bmp')))
    for i, fn in enumerate(files):
        print(i + 1)
        im = imamp * _im2double(_imread(fn))
        x = _put_slice(x, i, im)
    medim = np.mean(x, axis=2)

    imsub = medim - bgp
    imsub = imsub / np.max(imsub)          # normalize
    imsub1 = imsub * (imsub > 0)           # only above 0 data
    thr = np.mean(imsub1) + 2 * mstd(imsub1.ravel())  # threshold
    imsub2 = imsub * (imsub > thr)         # only above threshold
    # get approximate center
    fixed_results = None
    try:
        fixed_results = beadcentercalc(imsub2, 1, 0)
    except Exception as e:
        print(f'beadcentercalc failed for fixed trap: {e}')

    # %% get movable trap bead positions
    print('movable trap calc.. ')
    res = []
    results = None
    totimg = []
    for xi_ in xarr:
        for yj in yarr:
            files = sorted(glob.glob(os.path.join(imdir, f'*_{xi_}_{yj}_*.bmp')))
            print(f'doing...... X: {xi_} , Y: {yj}')
            # NOTE: as in the MATLAB code, the image stack `x` is re-used, so
            # only its first len(files) slices are overwritten here.
            for k, fn in enumerate(files):
                im = imamp * _im2double(_imread(fn))
                x = _put_slice(x, k, im)
            medim = np.mean(x, axis=2)
            imsub = medim - bgp                    # subtract background
            imsub = imsub / np.max(imsub)          # normalize
            imsub1 = imsub * (imsub > 0)           # only above 0 data
            thr = np.mean(imsub1) + 2 * mstd(imsub1.ravel())  # threshold
            imsub2 = imsub * (imsub > thr)         # only above threshold
            totimg.append(imsub2)
            # get approximate center
            try:
                results = beadcentercalc(imsub2, 0)
                res.append([float(xi_), float(yj), results.x0, results.y0,
                            results.pcorr, results.ssep, results.r2])
            except Exception:
                pass
    res = np.array(res, dtype=float)

    # %% keep positions within 2 std
    xpos = res[:, 2]
    idxx = np.where((xpos >= np.mean(xpos) - 2 * mstd(xpos)) & (xpos <= np.mean(xpos) + 2 * mstd(xpos)))[0]
    ypos = res[:, 3]
    idxy = np.where((ypos >= np.mean(ypos) - 2 * mstd(ypos)) & (ypos <= np.mean(ypos) + 2 * mstd(ypos)))[0]
    iid = np.intersect1d(idxx, idxy)

    fig, ax = plt.subplots()
    sc = ax.scatter(res[iid, 2], res[iid, 3], s=100, c=np.abs(res[iid, 5]))
    c = fig.colorbar(sc, ax=ax)
    c.set_label('sum squared error')
    ax.grid(True)
    ax.scatter(fixed_results.x0, fixed_results.y0, s=500, c='r', marker='x')
    ax.set_xlabel('X pixels')
    ax.set_ylabel('Y pixels')
    ax.set_xlim([np.min(res[iid, 2]) - 1, np.max(res[iid, 2]) + 15])
    ax.set_ylim([np.min(res[iid, 3]) - 1, np.max(res[iid, 3]) + 1])
    pts = np.column_stack([res[iid, 2], res[iid, 3]])
    fixed_x = float(griddata(pts, res[iid, 0], (fixed_results.x0, fixed_results.y0), method='cubic'))
    fixed_y = float(griddata(pts, res[iid, 1], (fixed_results.x0, fixed_results.y0), method='cubic'))
    print(fixed_x / 100)
    print(fixed_y / 100)
    print(np.sum(res[:, 5]))

    g1, = ax.plot(np.nan, np.nan, 'ko', label='movable trap pos. ')
    g2, = ax.plot(np.nan, np.nan, 'rx', label='fixed trap pos. ')
    ax.legend(handles=[g1, g2])
    ax.invert_yaxis()
    ax.set_aspect('equal', adjustable='datalim')
    ax.set_title(f'{strdate} {extra_title}\nfixed x(V): {fixed_x / 100:01.6f} fixed y(V): {fixed_y / 100:01.6f}')

    fig.savefig(os.path.join(imdir, imagesavename + '.png'))

    output = {}
    output['movable'] = res
    output[OUTPUT_FIXED_FIELD] = np.array([fixed_results.x0, fixed_results.y0])
    output['fixedV'] = np.array([fixed_x / 100, fixed_y / 100])

    # %% get beta that converts V to pixel. Y = bX formulation.
    order = np.lexsort((res[:, 1], res[:, 0]))  # sortrows(res, [1 2])
    sortedpos = res[order, :]

    VM = sortedpos[:, [0, 1]]
    VM1 = np.column_stack([VM / 100, np.ones(VM.shape[0])])  # include offset (values are V*100)
    PM = sortedpos[:, [2, 3]]

    # beta = PM'/VM1'  (least-squares solution of beta*VM1' = PM')
    beta = np.linalg.lstsq(VM1, PM, rcond=None)[0].T

    output['beta'] = beta

    # %% save results
    G.base_workspace['fixedxyres'] = output
    savemat(os.path.join(imdir, f'fixedxyres_{strdate}.mat'), output)

    # %%
    conditions = inputdlg(['Move Fixed images to new folder? (y or n)'], 'Move fixed x y images: ', ['y'])
    print(conditions)
    if conditions[0] == 'y':
        i = 0
        saveimgdir = os.path.join(imdir, f'fixedxyimdir_{i}')
        while os.path.isdir(saveimgdir):
            i = i + 1
            saveimgdir = os.path.join(imdir, f'fixedxyimdir_{i}')
        os.makedirs(saveimgdir)
        for fn in glob.glob(os.path.join(imdir, '*' + imagesavename + '.png')) + \
                glob.glob(os.path.join(imdir, '*' + imagesavename + '.fig')):
            shutil.move(fn, saveimgdir)
        for fn in glob.glob(os.path.join(imdir, '*.bmp')):
            shutil.move(fn, saveimgdir)
    else:
        print('did not move fixed xy images and figures. ')

    return output, res, results
