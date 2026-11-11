"""
Code_1  (Python version of Code_1.m)

XPD overview analysis: processivity / speed, burst statistics, lifetimes,
P_open and GC content plots.

Like the MATLAB file, this is a collection of cells ("# %%") that are meant
to be run one at a time (VS Code, Spyder, PyCharm or Jupyter all understand
"# %%" cells). Several cells re-define `outdtb` with different columns, so run
the data-loading cell that belongs to the plot you want first, just like in
MATLAB.

Column numbers of the `totdt` cell array are written 1-based (as in MATLAB)
through the helper `C(i, j)`; row index `i` is 0-based.

Dependencies that were not in the converted folder:
  * startup (defines analysisdir, processeddatadir, BoxPath) - set the paths
    in the "PATHS" block below, or keep a startup.py that defines them
  * MeanOccupancy_3 - only used in one cell
"""
# %% setup
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from scipy.optimize import curve_fit

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'New_Default_Analysis'))
from matlab_helpers import mstd, movmean, isempty  # noqa: E402

plt.close('all')

# ----------------------------------------------------------------- PATHS
try:
    from startup import analysisdir, processeddatadir, BoxPath  # your own startup.py
except ImportError:
    analysisdir = 'E:\\BOX\\Box\\UIUC_box_research\\progress_analysis_folder\\'   # <-- edit
    processeddatadir = 'E:\\BOX\\Box\\UIUC_box_research\\Data\\processeddata\\'   # <-- edit
    BoxPath = 'E:\\BOX\\Box'                                                    # <-- edit
# -----------------------------------------------------------------------

enddir = analysisdir + 'Recent\\overview\\'
os.makedirs(enddir, exist_ok=True)
figdir = enddir + 'workingfig\\'
os.makedirs(figdir, exist_ok=True)

TColr = np.array([
    [0.6667, 0.6667, 0.6667],
    [0.6667, 0.2667, 0.6000],
    [0.8000, 0.4000, 0.4667],
    [0.8667, 0.8000, 0.4667],
    [0.0667, 0.4667, 0.2000],
    [0.5333, 0.8000, 0.9333],
    [0.2000, 0.1333, 0.5333],
    [0.2667, 0.6667, 0.6000],
])

# long hairpin sequence used in the GC-content cells
SEQ = ('ATGCTAGCTGGGTGACCGCTGCCTAAAGCAAGTGAAGCAAGACCAGTGCCAATCGGTCCGAGCGGCCCTTTCAACGGTCATTGACCGACTTGTGTGGCAACTACTAGACGCGATGGACGACATGGACGCCGAAAAGCAGCATGACAAGGCCGTACAACAACCGCTAGGCTGCGTAGACGAGTGCCAGTTTCAATTGGTAGACACGCCGCTACAAAAAGTATCTAGGTGGGGCATTTAGGTCAGACACAAACAGTCCAGCTCAAAACCAAACGACCGACAGTGCGGACGGACAACGAACAATGCCAACTAAAGCTCAACCCAGGTGAATAGCGCCTCAAACCGGCCCGAGTCGTGCAATGGCCTGGTCTTCAACAGGACCGTACGGGTCGCCTTGTCGTGTGTGTACAGCGCCACCAATGCAGGCAGTGCAAGTGCGTAGTCCGCCTATAGCAATCGGGTGGGTCGTTTTAAGCCAAAAGACCGACTACCACGCTATCAGAAGTGGTACAGTTTGTAGGTGAGACGCCGCCAGTCCAGAAGACGACAGGGGGTGAACGACGGCGAGACTTAACGTCGTAGGCCAAAGTGGTGTCTTTCCAGCAAAAGACCGACCAGTCTCCTAAGCGGTCTTAAGAGACTGCTTATTAGAAAAGAAAAAAGAAAACATTATCACAGAAAACACAGGGGGACAAAACTCCCTATCGTTAGGGGGTTAAACTCCCTACAAAATAGGGAGCAAAATCCCCTAAAAGGGAGCAAAACTCCCTACGTGGTAAGACTCTACAAAAATAAACCAGGTTTGTACGGCGGAACGACGAACTATTATAAGTAAGACTGCTCAAGATTGAACCGAAGTAACGTGGCAAACTGTCCATTGAAACATTAGAGCGATTCAACTCTTAGCCACTAAGACAGGTAACCAAATAAGGTGGGTATCCAAAATGCGTCTTACCGTTCGTCGTGAAATTTGACAGCGAACCAGTCTAGACGCGGGCTTATTCGGAGTTCGTCGTATAAACTATCAGACCGCATTGGTAGTAGCTCTAGACGGTGTAATGCGAGGACAGGCCGTTTCAATGGAGACGGCTTCAACTCATAAAAACGACATAAACAGTATTACTGAGGACAACTATCTAGGTCATTACTGGAGTCTTGAGGTAGACCTAAACAAGTCTTGCGAGCCAACGGCGGCCCGCAAAAAATAACCACTCTTAGCGTCGTTGAACAGCGCGGTTAGCTCGGTACAGCAGCAGTTGCTGGGGGGTAAGTTCTTGTCGTTCGTCGTAACTCTTGAAACCTTAGGTCAGGGAGAAGGTGGACGACTAGACGCTGAATAGTTGCGGGTGTCGAAGGCGACAGAAGAGTCAAGGTTCGTAACGCTAAAACAATTCGTTGCGTGAGAGCTAAGCATCTCGGAGCAACGCAAACAAACGTGCTTGGTATACATTCATAAAGGAATCTATTGTTAACTAACTTACATACGTTTATTTACGTATGTGGTATCCACACCAAATTAAACTACGGGAAAAAGTCCCGACCTTACACATTCTCGCCCCAATAAATACGACAACAAAAAAACAATGAGCCCTTCCCGAAATGGAGAAGGCGTATTTGCGAAGGTAGTCGCAAATATCAATTTTTTTAGAAAGCCGGACGTACTTACCGGAACAACTAGCGCGAAACTATATGCGGCTCTAGAAATCGACAGAACCAAACGGG')


# =============================================================== helpers
def load_cell(path, varname):
    """load(path).varname for a MATLAB cell array -> numpy object array."""
    return loadmat(path, squeeze_me=False, struct_as_record=False)[varname]


def C(cellarr, i, j):
    """MATLAB cellarr{i+1, j}: row i is 0-based, column j is 1-based."""
    return cellarr[i, j - 1]


def num(x):
    """First element of a MATLAB value as float."""
    return float(np.asarray(x, dtype=float).ravel()[0])


def at(x, r, c):
    """MATLAB x(r, c) with 1-based r, c."""
    return float(np.atleast_2d(np.asarray(x, dtype=float))[r - 1, c - 1])


def lin(x, k):
    """MATLAB x(k) linear (column-major) 1-based index."""
    return float(np.asarray(x, dtype=float).ravel(order='F')[k - 1])


def cmean(x):
    return float(np.mean(np.asarray(x, dtype=float)))


def cmax(x):
    return float(np.max(np.asarray(x, dtype=float)))


def ksdensity(data, xi, bandwidth):
    """[f, xi] = ksdensity(data, xi, 'Bandwidth', bw) with a normal kernel."""
    data = np.asarray(data, dtype=float).ravel()
    xi = np.asarray(xi, dtype=float).ravel()
    z = (xi[:, None] - data[None, :]) / bandwidth
    f = np.exp(-0.5 * z**2).sum(axis=1) / (len(data) * bandwidth * np.sqrt(2 * np.pi))
    return f, xi


def ecdf_survivor(y):
    """[f, x] = ecdf(y, 'function', 'survivor')."""
    y = np.sort(np.asarray(y, dtype=float).ravel())
    y = y[~np.isnan(y)]
    xs, counts = np.unique(y, return_counts=True)
    F = np.cumsum(counts) / len(y)
    x = np.concatenate([[xs[0]], xs])
    f = np.concatenate([[1.0], 1 - F])
    return f, x


def gc_fraction(sub_seq):
    """oligoprop(sub_seq).GC / 100"""
    return (sub_seq.count('G') + sub_seq.count('C')) / len(sub_seq)


def new_fig(w=5.1 * 3, h=3.53 * 3):
    return plt.figure(figsize=(w, h))


def subplot_span(fig, m, n, idx):
    """MATLAB subplot(m, n, [idx...]) with 1-based, row-major indices."""
    idx = np.atleast_1d(idx) - 1
    rows, cols = idx // n, idx % n
    gs = GridSpec(m, n, figure=fig)
    return fig.add_subplot(gs[rows.min():rows.max() + 1, cols.min():cols.max() + 1])


def linkaxes(axs, which):
    for a in axs[1:]:
        if 'x' in which:
            a.sharex(axs[0])
        if 'y' in which:
            a.sharey(axs[0])


def style_ax(ax, fontsize=None, ygrid=True):
    """set(gca, 'Box','off', 'YGrid','on', 'XColor',[.3 .3 .3], 'color','none', 'LineWidth',3)"""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    for sp in ax.spines.values():
        sp.set_linewidth(3)
    ax.tick_params(width=3)
    ax.spines['bottom'].set_color((.3, .3, .3))
    ax.tick_params(axis='x', colors=(.3, .3, .3))
    if ygrid:
        ax.yaxis.grid(True)
    ax.set_facecolor('none')
    if fontsize:
        ax.tick_params(labelsize=fontsize)
        ax.xaxis.label.set_size(fontsize)
        ax.yaxis.label.set_size(fontsize)
        ax.title.set_size(fontsize)


def set_line_widths(ax, w=2):
    """set(findobj(gca,'type','line'),'linew',w)"""
    for ln in ax.lines:
        ln.set_linewidth(w)


def hide_xaxis(ax):
    """h.XAxis.Visible = 'off'"""
    ax.xaxis.set_visible(False)
    ax.spines['bottom'].set_visible(False)


def violin(ax, data, pos, color, bandwidth, show_data=True, alpha=0.3,
           width=0.3, markersize=24):
    """Violin(...) / one group of violinplot(...): KDE outline + jittered data."""
    data = np.asarray(data, dtype=float).ravel()
    data = data[~np.isnan(data)]
    if data.size == 0:
        return
    if data.size > 1 and np.ptp(data) > 0:
        yy = np.linspace(data.min(), data.max(), 200)
        d, _ = ksdensity(data, yy, bandwidth)
        d = d / d.max() * width
        ax.fill_betweenx(yy, pos - d, pos + d, facecolor=color, alpha=alpha,
                         edgecolor=color, linewidth=1)
    if show_data:
        if data.size > 1 and np.ptp(data) > 0:
            dd, _ = ksdensity(data, data, bandwidth)
            jit = (np.random.rand(data.size) - 0.5) * 2 * dd / dd.max() * width
        else:
            jit = np.zeros(data.size)
        ax.scatter(pos + jit, data, s=markersize, color=color, edgecolors='none')


def errorbar_msem(ax, msem, facecolor='w'):
    """errorbar(msem(:,1), msem(:,2), msem(:,3), 'ok', ... capsize 20, linewidth 3)"""
    return ax.errorbar(msem[:, 0], msem[:, 1], yerr=msem[:, 2], fmt='o', color='k',
                       markersize=10, markeredgecolor='k', markerfacecolor=facecolor,
                       capsize=10, elinewidth=3, capthick=3)


def group_by_site(tb, col, ldst, bin_by_value=True, maxval=None):
    """The repeated 'boxbin / boxdt / msem' loop.
    col is the 1-based MATLAB column of tb. If bin_by_value, boxbin holds the
    loading-site length (ldst(i)); otherwise the group index i."""
    boxbin, boxdt, scatterbin, avgarr, msem = [], [], [], [], []
    for i, L in enumerate(ldst, start=1):
        ttb = tb[tb[:, 0] == L, col - 1]
        if maxval is not None:
            ttb = ttb[ttb < maxval]
        if ttb.size > 0:
            boxbin += [L if bin_by_value else i] * len(ttb)
            boxdt += list(ttb)
            scatterbin += [i] * len(ttb)
            avgarr.append([i, np.mean(ttb)])
            msem.append([i, np.mean(ttb), mstd(ttb) / np.sqrt(len(ttb) + 1)])
        else:
            boxbin += [L if bin_by_value else i]
            boxdt += [0]
            scatterbin += [i]
            avgarr.append([i, np.nan])
            msem.append([i, np.nan, np.nan])
    return (np.array(boxbin, float), np.array(boxdt, float), np.array(scatterbin, float),
            np.array(avgarr, float), np.array(msem, float))


def load_popen():
    matfile = processeddatadir + 'Popen_5fkhp_SY_12pN.mat'
    p = loadmat(matfile, squeeze_me=True)

    class _P:
        pass
    popen = _P()
    popen.Hairpin_Position = np.asarray(p['Hairpin_Position'], dtype=float).ravel()
    popen.P_open = np.asarray(p['P_open'], dtype=float).ravel()
    return popen


def final_version_figure(boxdt, boxbin, msem, ldst, popen):
    """The 5x5 'Final version' figure (violins with broken y-axis, density and P_open)."""
    fig = new_fig()

    s1 = subplot_span(fig, 5, 5, [1, 2, 3])
    for iGroup, L in enumerate(ldst, start=1):
        violin(s1, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 10)
    errorbar_msem(s1, msem)
    s1.scatter(np.nan, np.nan, marker='o', color=(0.5, 0.5, 0.5))
    s1.set_title('1st burst processivity')
    s1.set_xticks(range(1, 8))
    s1.set_xticklabels([])
    s1.set_yticks(np.arange(150, 351, 50))
    s1.set_yticklabels(['', '200', '', '300', ''])
    s1.set_ylim([150, 350])
    set_line_widths(s1)
    style_ax(s1, 25)
    hide_xaxis(s1)

    s2 = subplot_span(fig, 5, 5, [6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23])
    for iGroup, L in enumerate(ldst, start=1):
        violin(s2, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 10)
    e1 = errorbar_msem(s2, msem)
    p1 = s2.scatter(np.nan, np.nan, marker='o', color=(0.5, 0.5, 0.5))
    s2.legend([p1, e1], ['data', r'mean $\pm$ sem'], loc='upper left')
    set_line_widths(s2)
    s2.set_ylabel(' Unwound(bp)')
    s2.set_xlabel('binding site length (nt)')
    s2.set_xticks(range(1, 8))
    s2.set_xticklabels(['10', '20', '30', '40', '50', '60', '70'])
    s2.set_ylim([0, 150])
    style_ax(s2, 25)

    s3 = subplot_span(fig, 5, 5, [4])
    xx = np.linspace(0, 350, 350)
    f, xi = ksdensity(boxdt, xx, 2)
    bins = np.arange(0, 351, 5)
    cts, edg = np.histogram(boxdt, bins)
    x = edg[:-1] + 0.5 * (edg[1] - edg[0])
    s3.barh(x, cts * np.max(f) / np.max(cts), height=edg[1] - edg[0],
            color=(0, 0, 0, 0.35), edgecolor='none')
    s3.plot(f, xi, 'k')
    s3.set_xticks(np.arange(0, 0.051, 0.025))
    s3.set_xticklabels([])
    s3.set_xlim([0, 0.05])
    s3.set_yticks(np.arange(150, 351, 50))
    s3.set_ylim([150, 350])
    s3.set_yticklabels([])
    s3.grid(True)
    set_line_widths(s3)
    style_ax(s3, 25)
    hide_xaxis(s3)

    s4 = subplot_span(fig, 5, 5, [9, 14, 19, 24])
    s4.barh(x, cts * np.max(f) / np.max(cts), height=edg[1] - edg[0],
            color=(0, 0, 0, 0.35), edgecolor='none')
    s4.plot(f, xi, 'k')
    s4.set_xticks(np.arange(0, 0.051, 0.025))
    s4.set_xticklabels(['0', '', '0.05'], rotation=0)
    s4.set_ylim([0, 150])
    s4.set_yticklabels([])
    s4.set_xlabel('Density')
    s4.grid(True)
    set_line_widths(s4)
    style_ax(s4, 25)

    s5 = subplot_span(fig, 5, 5, [5])
    s5.plot(popen.P_open, popen.Hairpin_Position, color=(1, 0.7, 0.7), linewidth=1, label='Popen')
    smooth_P_open = movmean(popen.P_open, 5)
    s5.plot(smooth_P_open, popen.Hairpin_Position, color=(1, 0, 0), linewidth=3, label='smooth 5')
    s5.set_yticks(np.arange(150, 351, 50))
    s5.set_ylim([150, 350])
    s5.set_xlim([0, 1])
    s5.grid(True)
    style_ax(s5, 25)
    s5.set_yticklabels([])
    s5.set_xticklabels([])
    hide_xaxis(s5)

    s6 = subplot_span(fig, 5, 5, [10, 15, 20, 25])
    s6.plot(popen.P_open, popen.Hairpin_Position, color=(1, 0.7, 0.7), linewidth=1, label='Popen')
    s6.plot(smooth_P_open, popen.Hairpin_Position, color=(1, 0, 0), linewidth=3, label='smooth 5')
    s6.set_ylim([0, 150])
    s6.set_xlim([0, 1])
    s6.grid(True)
    style_ax(s6, 25)
    s6.set_yticklabels([])
    s6.set_xlabel('P$_{open}$')

    linkaxes([s1, s2], 'x')
    linkaxes([s3, s4], 'x')
    linkaxes([s5, s6], 'x')
    linkaxes([s1, s3, s5], 'y')
    linkaxes([s2, s4, s6], 'y')
    return fig


def gc_profile(seq, seq_len=10, with_pupy=False):
    GC = np.full(len(seq), np.nan)
    PyPu = np.full((2, len(seq)), np.nan)
    for i in range(len(seq)):
        st = i
        ed = min(len(seq) - 1, i + seq_len)
        sub_seq = seq[st:ed + 1]
        GC[i] = gc_fraction(sub_seq)
        if with_pupy:
            pu = sub_seq.count('A') + sub_seq.count('G')
            py = sub_seq.count('C') + sub_seq.count('T')
            total_bases = len(sub_seq)
            PyPu[0, i] = pu / total_bases  # purine fraction
            PyPu[1, i] = py / total_bases  # pyrimidine fraction
    return GC, PyPu


def gc_figure(popen, GC, PyPu=None):
    P_open = popen.P_open
    x = popen.Hairpin_Position
    n = len(x)
    fig = plt.figure()
    s1 = subplot_span(fig, 1, 4, [1, 2, 3])
    s1.plot(x + 1, GC[:n], color=(0.7, 0.7, 1))
    s1.plot(x, P_open, color=(1, 0.7, 0.7))
    s1.plot(x + 1, movmean(GC[:n], 10), linewidth=3, color='b', label='%GC')
    s1.plot(x, movmean(P_open, 10), linewidth=3, color='r', label='P$_{open}$')
    if PyPu is not None:
        s1.plot(x, movmean(PyPu[0, :n], 10), linewidth=3, color=(0.1, 1, 0.1), label='pur')
        s1.plot(x, movmean(PyPu[1, :n], 10), linewidth=3, color=(1, 0, 1), label='pyr')
    s1.set_xlim([1, 75])
    s1.set_ylim([0, 1])
    s1.grid(True)
    s1.set_xlabel('Hairpin Position(bp)')
    s1.set_ylabel('%GC & P$_{open}$')
    s1.legend(loc='upper left')
    style_ax(s1, 25)

    s2 = subplot_span(fig, 1, 4, [4])
    s2.plot(x, GC[:n], color=(0.8, 0.8, 1))
    s2.plot(x, P_open, color=(1, 0.7, 0.7))
    s2.plot(x, movmean(GC[:n], 10), linewidth=3, color='b')
    s2.plot(x, movmean(P_open, 10), linewidth=3, color='r')
    s2.set_xlim([75, 350])
    s2.set_ylim([0, 1])
    linkaxes([s1, s2], 'y')
    s2.set_yticklabels([])
    s2.grid(True)
    style_ax(s2, 25, ygrid=False)
    return fig


# %% load data (line 30)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
fhi = 13
flo = 11.5

# %% outdtb: [loading site, mean force, end of activity, max extension, 1st burst] (line 35)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
fhi = 13
flo = 11.5
outdtb = []
for i in range(totdt.shape[0]):
    if not isempty(C(totdt, i, 16)) and flo <= cmean(C(totdt, i, 14)) <= fhi:
        if not isempty(C(totdt, i, 8)):
            outdtb.append([num(C(totdt, i, 4)), cmean(C(totdt, i, 14)), num(C(totdt, i, 16)),
                           cmax(C(totdt, i, 7)), at(C(totdt, i, 8), 1, 2)])
outdtb = np.array(outdtb)

# %% outdtb: [loading site, processivity, speed] (line 63)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
fhi = 13
flo = 11.5
outdtb = []
for i in range(totdt.shape[0]):
    if not isempty(C(totdt, i, 16)) and flo <= cmean(C(totdt, i, 14)) <= fhi:
        if not isempty(C(totdt, i, 8)):
            outdtb.append([num(C(totdt, i, 4)), at(C(totdt, i, 8), 1, 2),
                           at(C(totdt, i, 8), 1, 2) / at(C(totdt, i, 8), 1, 1)])
outdtb = np.array(outdtb)

# %% processivity vs speed
dt = outdtb
plt.figure()
plt.plot(dt[:, 1], dt[:, 2], '.', color=(0.5, 0.5, 0.5), markersize=35)
plt.xlabel(' Processivity')
plt.ylabel('Brust Unwinding Speed(bp/s)')
plt.grid(True)
figname = 'proc_v_speed_1'
# plt.savefig(figdir + figname + '.png')

# %% hist3
dtt = np.column_stack([dt[:, 1], dt[:, 2]])
ctrs_x = np.arange(0, 301, 10)
ctrs_y = np.arange(0, 26, 1)
ex_ = np.concatenate([[-np.inf], (ctrs_x[:-1] + ctrs_x[1:]) / 2, [np.inf]])
ey_ = np.concatenate([[-np.inf], (ctrs_y[:-1] + ctrs_y[1:]) / 2, [np.inf]])
H, _, _ = np.histogram2d(dtt[:, 0], dtt[:, 1], bins=[ex_, ey_])
plt.figure()
plt.pcolormesh(ctrs_x, ctrs_y, H.T, shading='nearest')
plt.xlabel('MPG')
plt.ylabel('Weight')
plt.colorbar()

# %% binned mean / std / ste
x = dt[:, 1]
y = dt[:, 2]
idx = np.argsort(x, kind='stable')
sx = x[idx]
sy = y[idx]

l = 20
n = int(np.ceil(len(sx) / l))
errarr = []
for i in range(1, n + 1):
    stidx = (i - 1) * l + 1
    edidx = min(i * l, len(sx))
    cx = sx[stidx - 1:edidx]
    cy = sy[stidx - 1:edidx]
    mx = np.mean(cx); stx = mstd(cx); stex = mstd(cx) / len(cx)
    my = np.mean(cy); sty = mstd(cy); stey = mstd(cy) / len(cy)
    errarr.append([mx, stx, stex, my, sty, stey])
errarr = np.array(errarr)

# %%
plt.figure()
ex = errarr[:, 0]; ey = errarr[:, 3]
esx = errarr[:, 2]; esy = errarr[:, 5]  # ste
plt.errorbar(ex, ey, yerr=esy, xerr=esx, fmt='o', color='k', elinewidth=3, capsize=10,
             markersize=10, markeredgecolor='r')
plt.xlabel('Processivity(bp)')
plt.ylabel('Burst Unwinding Speed(bp/s)')
plt.grid(True)
plt.title('Proc and Speed dependence (STE) ')
figname = 'proc_v_speed_ste'

# %% mean occupancy (needs MeanOccupancy_3)
from MeanOccupancy_3 import MeanOccupancy_3  # noqa: E402  (not part of the converted folder)
meannoXPD = []
ldst = np.arange(10, 71, 1)
a = 10
probout = np.zeros((6, 7))
for L in ldst:
    E = -np.log(5 * 10 / L)
    Nav, Cc, prob = MeanOccupancy_3(L, a, E)
    meannoXPD.append([L, Nav])
meannoXPD = np.array(meannoXPD)

# %%
plt.figure()
plt.plot(meannoXPD[:, 0], meannoXPD[:, 1], 'bo--')
plt.grid(True)

# %% P_open
popen = load_popen()
plt.figure()
plt.plot(popen.Hairpin_Position, popen.P_open, 'b', label='Popen')
smooth_P_open = movmean(popen.P_open, 5)
plt.plot(popen.Hairpin_Position, smooth_P_open, 'r--', label='smooth 5')

# %% P_open histograms
plt.figure()
ax = plt.gca()
bw = 0.05
bins = np.arange(0, 1 + bw / 2, bw)

dtp = popen.P_open
cts, edg = np.histogram(dtp, bins)
x = edg[:-1] + 0.5 * (edg[1] - edg[0])
ax.bar(x, cts / np.sum(cts), width=bw, facecolor=(0, 0, 0, 0.4), edgecolor='k', linewidth=4, label='All')

dtp = popen.P_open[:26]
cts, edg = np.histogram(dtp, bins)
x = edg[:-1] + 0.5 * (edg[1] - edg[0])
ax.bar(x, cts / np.sum(cts), width=bw, facecolor=(1, 0, 0, 0.4), edgecolor=(0.8, 0.1, 0.1), linewidth=4,
       label='1st barrier')
ax.legend()
ax.set_ylabel('Density')
ax.set_xlabel('Loading Site Length (nt)')

# %% 1st burst processivity, violinplot + histogram (line 220)
tb = outdtb
ldst = np.arange(10, 71, 10)
boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 5, ldst, bin_by_value=False)
fig = new_fig()
ax = subplot_span(fig, 1, 4, [1, 2, 3])
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for k, g in enumerate(np.unique(boxbin)):
    violin(ax, boxdt[boxbin == g], k + 1, colors[k % len(colors)], 10)
errorbar_msem(ax, msem)
ax.legend(['data', r'mean $\pm$ sem'], loc='upper left')
set_line_widths(ax)
ax.set_ylabel(' Unwound(bp)')
ax.set_xlabel('binding site length (nt)')
ax.set_title('1st burst processivity')
ax.set_xticks(range(1, 8))
ax.set_xticklabels(['10', '20', '30', '40', '50', '60', '70'])
ax.set_ylim([0, 350])
style_ax(ax, 20)

ax = subplot_span(fig, 1, 4, [4])
bins = np.concatenate([np.arange(0, 351, 5), [np.inf]])
cts, edg = np.histogram(boxdt, bins)
x = edg[:-1] + 0.5 * (edg[1] - edg[0])
ax.barh(x, cts, height=edg[1] - edg[0], facecolor=(0, 0, 0, 0.1), edgecolor='k', linewidth=4)
ax.plot(popen.P_open * np.max(cts) / np.max(popen.P_open), popen.Hairpin_Position, 'r', label='Popen')
smooth_P_open = movmean(popen.P_open, 5)
ax.plot(smooth_P_open * np.max(cts) / np.max(smooth_P_open), popen.Hairpin_Position, 'b', label='smooth 5')
ax.set_ylim([0, 350])
ax.grid(True)
set_line_widths(ax)
style_ax(ax, 20)
figname = 'burst1unwound_with_hist'

# %% 1st burst processivity, Violin + histogram + P_open on top axis (line 314)
tb = outdtb
ldst = np.arange(10, 71, 10)
boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 5, ldst)
fig = new_fig()
ax = subplot_span(fig, 1, 4, [1, 2, 3])
for iGroup, L in enumerate(ldst, start=1):
    violin(ax, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 10)
e1 = errorbar_msem(ax, msem)
p1 = ax.scatter(np.nan, np.nan, marker='o', color=(0.5, 0.5, 0.5))
ax.legend([p1, e1], ['data', r'mean $\pm$ sem'], loc='upper left')
set_line_widths(ax)
ax.set_ylabel(' Unwound(bp)')
ax.set_xlabel('binding site length (nt)')
ax.set_title('1st burst processivity')
ax.set_xticks(range(1, 8))
ax.set_xticklabels(['10', '20', '30', '40', '50', '60', '70'])
ax.set_ylim([0, 350])
style_ax(ax, 25)

ax = subplot_span(fig, 1, 4, [4])
bins = np.arange(0, 351, 5)
cts, edg = np.histogram(boxdt, bins)
x = edg[:-1] + 0.5 * (edg[1] - edg[0])
ax.barh(x, cts, height=edg[1] - edg[0], facecolor=(0, 0, 0, 0.1), edgecolor='k', linewidth=4)
ax.set_xticks(np.arange(0, 81, 20))
ax.set_xlim([0, 75])
ax.set_ylim([0, 350])
ax.set_xlabel('counts')
hax2 = ax.twiny()  # second x axis on top for P_open
hax2.plot(popen.P_open, popen.Hairpin_Position, 'r', label='Popen')
smooth_P_open = movmean(popen.P_open, 5)
hax2.plot(smooth_P_open, popen.Hairpin_Position, 'b', label='smooth 5')
hax2.set_xlim([0, 1])
hax2.set_xlabel('P$_{open}$')
ax.set_ylim([0, 350])
ax.grid(True)
set_line_widths(ax)
set_line_widths(hax2)
style_ax(ax, 25)
hax2.tick_params(labelsize=25)
hax2.xaxis.label.set_size(25)
figname = 'burst1unwound_with_hist'

# %% separate out things (line 456)
tb = outdtb
ldst = np.arange(10, 71, 10)
boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 5, ldst)
fig = new_fig()

s1 = subplot_span(fig, 5, 5, [1, 2, 3])
for iGroup, L in enumerate(ldst, start=1):
    violin(s1, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 10)
errorbar_msem(s1, msem)
s1.set_title('1st burst processivity')
s1.set_xticks(range(1, 8))
s1.set_xticklabels([])
s1.set_yticks(np.arange(150, 351, 50))
s1.set_yticklabels(['', '200', '', '300', ''])
s1.set_ylim([150, 350])
set_line_widths(s1)
style_ax(s1, 25)
hide_xaxis(s1)

s2 = subplot_span(fig, 5, 5, [6, 7, 8, 11, 12, 13, 16, 17, 18, 21, 22, 23])
for iGroup, L in enumerate(ldst, start=1):
    violin(s2, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 10)
e1 = errorbar_msem(s2, msem)
p1 = s2.scatter(np.nan, np.nan, marker='o', color=(0.5, 0.5, 0.5))
s2.legend([p1, e1], ['data', r'mean $\pm$ sem'], loc='upper left')
set_line_widths(s2)
s2.set_ylabel(' Unwound(bp)')
s2.set_xlabel('binding site length (nt)')
s2.set_xticks(range(1, 8))
s2.set_xticklabels(['10', '20', '30', '40', '50', '60', '70'])
s2.set_ylim([0, 150])
style_ax(s2, 25)

s3 = subplot_span(fig, 5, 5, [4])
bins = np.arange(0, 351, 5)
cts, edg = np.histogram(boxdt, bins)
x = edg[:-1] + 0.5 * (edg[1] - edg[0])
s3.barh(x, cts, height=edg[1] - edg[0], facecolor=(0, 0, 0, 0.1), edgecolor='k', linewidth=4)
s3.set_xticks(np.arange(0, 81, 30))
s3.set_xticklabels([])
s3.set_yticks(np.arange(150, 351, 50))
s3.set_ylim([150, 350])
s3.set_yticklabels([])
s3.set_xlim([0, 0.1])  # (as in the original, this overrides the earlier xlim([0 75]))
s3.grid(True)
set_line_widths(s3)
style_ax(s3, 25)
hide_xaxis(s3)

s4 = subplot_span(fig, 5, 5, [9, 14, 19, 24])
s4.barh(x, cts / np.sum(cts), height=edg[1] - edg[0], facecolor=(0, 0, 0, 0.1), edgecolor='k', linewidth=4)
xx = np.linspace(0, 350, 350)
f, xi = ksdensity(boxdt, xx, 3)
s4.plot(f, xi, 'k')
s4.set_ylim([0, 150])
s4.set_xlabel('counts')
s4.grid(True)
set_line_widths(s4)
style_ax(s4, 25)

s5 = subplot_span(fig, 5, 5, [5])
s5.plot(popen.P_open, popen.Hairpin_Position, color=(1, 0.7, 0.7), linewidth=1, label='Popen')
smooth_P_open = movmean(popen.P_open, 5)
s5.plot(smooth_P_open, popen.Hairpin_Position, color=(1, 0, 0), linewidth=3, label='smooth 5')
s5.set_yticks(np.arange(150, 351, 50))
s5.set_ylim([150, 350])
s5.set_xlim([0, 1])
s5.grid(True)
style_ax(s5, 25)
s5.set_yticklabels([])
s5.set_xticklabels([])
hide_xaxis(s5)

s6 = subplot_span(fig, 5, 5, [10, 15, 20, 25])
s6.plot(popen.P_open, popen.Hairpin_Position, color=(1, 0.7, 0.7), linewidth=1, label='Popen')
s6.plot(smooth_P_open, popen.Hairpin_Position, color=(1, 0, 0), linewidth=3, label='smooth 5')
s6.set_ylim([0, 150])
s6.set_xlim([0, 1])
s6.grid(True)
style_ax(s6, 25)
s6.set_yticklabels([])
s6.set_xlabel('P$_{open}$')

linkaxes([s1, s2], 'x')
linkaxes([s3, s4], 'x')
linkaxes([s5, s6], 'x')
linkaxes([s1, s3, s5], 'y')
linkaxes([s2, s4, s6], 'y')
figname = 'burst1unwound_with_hist'

# %% Final version (line 790)
tb = outdtb
ldst = np.arange(10, 71, 10)
boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 5, ldst)
final_version_figure(boxdt, boxbin, msem, ldst, popen)
figname = 'burst1unwound_with_hist'

# %% density as area plot (line 1138)
bins = np.arange(0, 351, 1)
cts, edg = np.histogram(boxdt, bins)
x = edg[:-1] + 0.5 * (edg[1] - edg[0])
xx = np.linspace(0, 350, 350)
f, xi = ksdensity(boxdt, xx, 3)
plt.figure()
plt.fill_between(xi, f, 1)  # area(xi, f, 1)

# %% max # XPD vs ssDNA (line 1171)
MParr = np.array([[30, 2], [40, 3], [50, 4], [60, 4], [70, 5], [80, 5]])
SYMP = np.array([[10, 1], [20, 1], [30, 2], [40, 3], [50, 3], [60, 4], [70, 5], [80, 6]])

plt.figure()
x = MParr[:, 0]; y = MParr[:, 1]
p = np.polyfit(x, y, 1); print(p)
xx = np.linspace(10, 80, 100)
yy = np.polyval(p, xx)
plt.plot(x, y, 'ro')
plt.plot(xx, yy, 'r-')
x = SYMP[:, 0]; y = SYMP[:, 1]
p = np.polyfit(x, y, 1); print(p)
yy = np.polyval(p, xx)
plt.plot(SYMP[:, 0], SYMP[:, 1], 'bx')
plt.plot(xx, yy, 'b-')
plt.grid(True)
plt.ylim([0, 7])

plt.figure()
y = MParr[:, 0]; x = MParr[:, 1]
p = np.polyfit(x, y, 1); print(p)
xx = np.linspace(1, 7, 100)
yy = np.polyval(p, xx)
plt.plot(x, y, 'ro')
plt.plot(xx, yy, 'r-')
y = SYMP[:, 0]; x = SYMP[:, 1]
p = np.polyfit(x, y, 1); print(p)
yy = np.polyval(p, xx)
plt.plot(x, y, 'bx')
plt.plot(xx, yy, 'b-')
plt.grid(True)
plt.ylim([0, 90])
plt.xlim([0, 7])
plt.xlabel('MAX # XPD')
plt.ylabel('ssDNA')

# %% (line 1220)
plt.figure()
x = MParr[:, 0]; y = MParr[:, 1]
p = np.polyfit(x, y, 1); print(p)
xx = np.linspace(10, 80, 100)
yy = np.polyval(p, xx)
plt.plot(x, y, 'ro')
plt.plot(xx, yy, 'r-')
x = SYMP[:, 0]; y = SYMP[:, 1]
p = np.polyfit(x, y, 1); print(p)
yy = np.polyval(p, xx)
plt.plot(SYMP[:, 0], SYMP[:, 1], 'bx')
plt.plot(xx, yy, 'b-')
plt.grid(True)
plt.ylim([0, 7])

# %% (line 1248)
plt.figure()
y = MParr[:, 0]; x = MParr[:, 1]
p = np.polyfit(x, y, 1); print(p)
xx = np.linspace(0, 7, 100)
yy = np.polyval(p, xx)
plt.plot(x, y, 'd', markersize=35, markerfacecolor='r', markeredgecolor='none', label='Spies Group values')
plt.plot(xx, yy, 'r-', label=f'slope = {p[0]:02.1f}')
y = SYMP[:, 0]; x = SYMP[:, 1]
p = np.polyfit(x, y, 1); print(p)
yy = np.polyval(p, xx)
plt.plot(x, y, 's', markersize=35, markerfacecolor='b', markeredgecolor='none', label='My values')
plt.plot(xx, yy, 'b-', label=f'slope = {p[0]:02.1f}')
plt.grid(True)
plt.ylim([0, 90])
plt.xlim([0, 7])
plt.xlabel('MAX # XPD')
plt.ylabel('ssDNA')
plt.legend()

# %% GC content (line 1277)
GC, _ = gc_profile(SEQ, 10)
gc_figure(popen, GC)

# %% GC content (line 1339)
GC, _ = gc_profile(SEQ, 10)
gc_figure(popen, GC)

# %% GC content + purine / pyrimidine (line 1404)
GC, PyPu = gc_profile(SEQ, 10, with_pupy=True)
gc_figure(popen, GC, PyPu)

# %% outdtb with all bursts (line 1489)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
totdt = load_cell(processeddatadir + 'totdt_endofactivity_210503.mat', 'out')
fhi = 13
flo = 11.5
outdtb = []
for i in range(totdt.shape[0]):
    if not isempty(C(totdt, i, 16)) and flo <= cmean(C(totdt, i, 14)) <= fhi:
        if not isempty(C(totdt, i, 8)):
            b = np.atleast_2d(np.asarray(C(totdt, i, 8), dtype=float))[:, 1]
            outdtb.append(np.concatenate([[num(C(totdt, i, 4)), cmean(C(totdt, i, 14))],
                                          b, np.full(100 - len(b), np.nan)]))
outdtb = np.array(outdtb)

# %% distribution of burst sizes per burst number (line 1516)
ldstarr = np.array([10, 20, 30, 40, 50, 70])
df = outdtb
ldst = np.array([10, 40, 20, 50, 30, 70])
fig = plt.figure(figsize=(5.1 * 2.5, 3.53 * 2.5))
for i in range(len(ldst)):
    ax = fig.add_subplot(3, 2, i + 1)
    dff = df[df[:, 0] == ldst[i], 2:7]
    boxbin = []
    boxdt = []
    msem = []
    for ii in range(1, 6):
        dfff = dff[:, ii - 1]
        dfff = dfff[~np.isnan(dfff)]
        boxbin += [ii] * len(dfff)
        boxdt += list(dfff)
        msem.append([ii, np.mean(dfff) if len(dfff) else np.nan,
                     mstd(dfff) / np.sqrt(len(dfff)) if len(dfff) else np.nan])
    boxbin = np.array(boxbin)
    boxdt = np.array(boxdt)
    msem = np.array(msem)
    col = TColr[int(ldst[i] / 10) - 1]
    for k, g in enumerate(np.unique(boxbin)):
        violin(ax, boxdt[boxbin == g], k + 1, col, 10, alpha=0.3, markersize=15)
    vs2, = ax.plot(np.nan, np.nan, '.', markersize=15, color=col)
    er = ax.errorbar(msem[:, 0], msem[:, 1], yerr=msem[:, 2], fmt='x', color='k', markersize=10,
                     capsize=10, elinewidth=2)
    errorbar_msem(ax, msem)
    # (in MATLAB the second legend() call replaces the first one)
    ax.legend([vs2], [f'{ldst[i]}dt'], loc='upper right')
    if i == 2:
        ax.set_ylabel('Bp unwound')
    if i >= 4:
        ax.set_xlabel('Burst #')
    ax.set_xlim([0, 6])
    ax.set_ylim([0, 350])
    ax.set_yticks(np.arange(0, 351, 50))
    ax.set_yticklabels(['0', '', '100', '', '200', '', '300', ''])
    ax.grid(True)
figname = 'dist_burstgen_1'

# %% Distribution of bursts (line 1598)
burstarr = []
for i in range(len(ldst)):
    dff = df[df[:, 0] == ldst[i], :]
    numarr = []
    for ii in range(dff.shape[0]):
        dfff = dff[ii, :]
        curnum = np.sum(~np.isnan(dfff))
        numarr.append(curnum)
    numarr = np.array(numarr, dtype=float)
    m = np.mean(numarr)
    sem = mstd(numarr) / np.sqrt(len(numarr))
    burstarr.append([ldst[i], m, sem])
burstarr = np.array(burstarr)


# %% lifetime violin (line 1626)  -- needs the outdtb of the line-35 cell
def lifetime_violin(outdtb, facecolor):
    tb = outdtb
    ldst = np.arange(10, 71, 10)
    boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 3, ldst)
    idx = np.where(boxdt < 200)[0]
    boxdt = boxdt[idx]
    boxbin = boxbin[idx]
    fig = plt.figure(figsize=(5.1 * 2.5, 3.53 * 2.5))
    ax = fig.gca()
    for iGroup, L in enumerate(ldst, start=1):
        violin(ax, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 15, show_data=False)
    errorbar_msem(ax, msem, facecolor=facecolor)
    set_line_widths(ax)
    ax.set_ylabel(' Activity duration(s) ')
    ax.set_xlabel('loading site length(nt)')
    ax.set_title('Lifetime')
    ax.set_xticks(range(1, 8))
    ax.set_xticklabels(['10', '20', '30', '40', '50', '60', '70'])
    ax.set_ylim([0, 200])
    style_ax(ax, 35)
    return fig


lifetime_violin(outdtb, 'k')
figname = 'lifetime_violin_2'

# %% outdtb for lifetimes (line 1725)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
fhi = 140
flo = 0.11
outdtb = []
for i in range(totdt.shape[0]):
    if not isempty(C(totdt, i, 16)) and flo <= cmean(C(totdt, i, 14)) <= fhi:
        if not isempty(C(totdt, i, 8)):
            burstmax = num(C(totdt, i, 16))
            outdtb.append([num(C(totdt, i, 4)), cmean(C(totdt, i, 14)), num(C(totdt, i, 16)),
                           cmax(C(totdt, i, 7)), at(C(totdt, i, 8), 1, 2),
                           min(num(C(totdt, i, 16)), lin(C(totdt, i, 9), 2))])
outdtb = np.array(outdtb)

# %% survival curves (line 1755)
tb = outdtb
ldst = np.arange(10, 71, 10)
plt.figure()
ax = plt.gca()
for i in range(len(ldst)):
    ttb = tb[tb[:, 0] == ldst[i], 2]
    if ttb.size > 0:
        ttb = ttb[ttb < 200]
        f, x = ecdf_survivor(ttb)
        ax.plot(x, f, '-', color=TColr[i], linewidth=3, markersize=20, label=f'{ldst[i]}-dT')
ax.legend()

ttb = tb[tb[:, 0] == 10, 2]
f, x = ecdf_survivor(ttb)
ok = np.isfinite(x) & np.isfinite(f)  # prepareCurveData
xData, yData = x[ok], f[ok]
ft = lambda x, a, b: a * np.exp(-b * x)  # 'a*exp(-b*(x))'
popt, _ = curve_fit(ft, xData, yData, p0=[0.6668284765108, 0.0170662910249856])
xl = ax.get_xlim()
xf = np.linspace(xl[0], xl[1], 200)
ax.plot(xf, ft(xf, *popt), 'r',
        label=f'koff: {popt[1]:03.3f}' + r'$\pm$ 0.0005 s$^{-1}$')
ax.set_ylim([0, 1])
ax.legend()
ax.set_xlabel(' activity duration(s) ')
ax.set_ylabel('survival prob')
ax.set_title('Lifetime')
style_ax(ax, 35)
figname = 'lifetime_survival'

# %% lifetime violin with ttb < 200 applied per group (line 1814)
tb = outdtb
ldst = np.arange(10, 71, 10)
boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 3, ldst, maxval=200)
idx = np.where(boxdt < 200)[0]
boxdt = boxdt[idx]
boxbin = boxbin[idx]
fig = plt.figure(figsize=(5.1 * 2.5, 3.53 * 2.5))
ax = fig.gca()
for iGroup, L in enumerate(ldst, start=1):
    violin(ax, boxdt[boxbin == L], iGroup, TColr[iGroup - 1], 15, show_data=False)
errorbar_msem(ax, msem, facecolor=(1, 1, 1))
set_line_widths(ax)
ax.set_ylabel(' Activity duration(s) ')
ax.set_xlabel('loading site length(nt)')
ax.set_title('Lifetime')
ax.set_xticks(range(1, 8))
ax.set_xticklabels(['10', '20', '30', '40', '50', '60', '70'])
ax.set_ylim([0, 200])
style_ax(ax, 35)
figname = 'lifetime_violin_2'

# %% strand switch / backslide data (line 1915)
processeddir = 'E:\\Box\\Box\\UIUC_box_research\\progress_analysis_folder\\Past\\meeting200826\\'
totdt = load_cell(processeddir + 'out_200819_strandswitchorbackslide.mat', 'out')
out = totdt
outdt = []
for i in range(out.shape[0]):
    if not isempty(C(out, i, 17)):
        if lin(C(out, i, 17), 1) == 1:
            outdt.append([num(C(out, i, 4)), lin(C(out, i, 17), 2), lin(C(out, i, 17), 3)])
outdt = np.array(outdt)

# %% Final version (strand switch) (line 1933)
tb = outdt
ldst = np.arange(10, 71, 10)
boxbin, boxdt, scatterbin, avgarr, msem = group_by_site(tb, 3, ldst)
final_version_figure(boxdt, boxbin, msem, ldst, popen)

# %% (line 2279)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
fhi = 13
flo = 11.5
outdtb = []
for i in range(totdt.shape[0]):
    if not isempty(C(totdt, i, 16)) and flo <= cmean(C(totdt, i, 14)) <= fhi:
        if not isempty(C(totdt, i, 8)):
            outdtb.append([num(C(totdt, i, 4)), cmean(C(totdt, i, 14)), num(C(totdt, i, 16)),
                           cmax(C(totdt, i, 7)), at(C(totdt, i, 8), 1, 2)])
outdtb = np.array(outdtb)

# %% (line 2309)
processeddir = BoxPath + '\\UIUC_box_research\\progress_analysis_folder\\Past\\meeting200826\\'
totdt = load_cell(processeddir + 'out_200819_strandswitchorbackslide.mat', 'out')
out = totdt
outdt = []
fhi = 13
flo = 11.5
for i in range(out.shape[0]):
    if (not isempty(C(out, i, 17)) and not isempty(C(totdt, i, 16))
            and flo <= cmean(C(totdt, i, 14)) <= fhi):
        if lin(C(out, i, 17), 1) == 1:
            outdt.append([num(C(out, i, 4)), lin(C(out, i, 17), 2), lin(C(out, i, 17), 3),
                          at(C(totdt, i, 8), 1, 2)])
outdt = np.array(outdt)

# %% (line 2328)
plt.figure()
tb = outdt
ldst = np.arange(10, 71, 10)
for i in range(len(ldst)):
    ttb = tb[tb[:, 0] == ldst[i]][:, [3, 2]]
    plt.scatter(ttb[:, 0], ttb[:, 1], s=50, facecolors='none', edgecolors=[TColr[i]])

# %% (line 2357)
plt.figure()
tb = outdt
for i in range(len(ldst)):
    ttb = tb[tb[:, 0] == ldst[i]][:, [3, 2]]
    plt.scatter(ttb[:, 0], ttb[:, 1], s=50, facecolors='none', edgecolors=[TColr[i]])

# %% first strand switching position vs processivity (line 2387)
totdt = load_cell(processeddatadir + 'totdt_endofactivity_200331_2.mat', 'totdt')
fhi = 13
flo = 11.5
outdtb = []
for i in range(totdt.shape[0]):
    if not isempty(C(totdt, i, 16)) and flo <= cmean(C(totdt, i, 14)) <= fhi:
        if not isempty(C(totdt, i, 8)):
            outdtb.append([num(C(totdt, i, 4)), cmean(C(totdt, i, 14)), num(C(totdt, i, 16)),
                           cmax(C(totdt, i, 7)), at(C(totdt, i, 8), 1, 2), at(C(totdt, i, 10), 1, 2)])
outdtb = np.array(outdtb)

# %% (line 2415)
dtarr = np.arange(10, 71, 10)
plt.figure()
for i in range(len(dtarr)):
    ndt = dtarr[i]
    tdt = outdtb[outdtb[:, 0] == ndt][:, [5, 4]]
    plt.plot(tdt[:, 0], tdt[:, 1], 'o', markersize=10, color=TColr[i], markerfacecolor='none',
             label=f'{ndt} dT')
plt.xlabel('Frist Strand switching position (bp)')
plt.ylabel('Processivity (bp)')
plt.xlim([3, 400])
plt.ylim([3, 400])
plt.legend()

plt.show()
