"""
d240501_XPD_unwinding_data_generation_for_bnp_step

Simulate XPD unwinding steps (Zhi's kinetic model), make an evenly sampled
trace, add white Gaussian noise and save a downsampled noisy trace.

Needs ModelHairpinPopen and SimKineticCycleSequenceTest (Zhi's simulation
code, not part of the converted folder) to be importable as Python modules.
Run cell by cell (# %%) in VS Code / Spyder, or as a whole script.
"""
# %%
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'New_Default_Analysis'))
from matlab_helpers import downsample  # noqa: E402

rng = np.random.RandomState(1231)  # rng(1231, "twister")


def stairs_xy(x, y):
    """[xb, yb] = stairs(x, y) for vectors."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    xb = np.repeat(x, 2)[1:]
    yb = np.repeat(y, 2)[:-1]
    return xb, yb


def awgn(sig, snr_db):
    """awgn(sig, snr): add white Gaussian noise assuming 0 dBW signal power."""
    sig = np.asarray(sig, dtype=float)
    noise_power = 10 ** (-snr_db / 10)
    return sig + np.sqrt(noise_power) * rng.standard_normal(sig.shape)


# %% Input
os.chdir(r'D:\Box_D\Box\UIUC_box_research\Code\ResearchProject\code_from_other_people\ZhiCode\kF_kB_Model (Zhis simulation codes)')
sys.path.insert(0, os.getcwd())
Path_Results_1 = r'D:\Box_D\Box\UIUC_box_research\progress_analysis_folder\Recent\meeting240501_zhisim_bnp_step/'

Seq = 1
ATP = 500           # uM
F_feedback = 12     # pN
NUM_Simulation = 100

# simulation parameters
# k is a list of length N where N is the number of states. Each k[i] is an
# array of the transition rates OUT of that state i
# k{1} = [keff_on*ATP k_0]; k{2} = [k_add k_inf k_minus5]; k{3} = [k_add5]; % new
k = [np.array([0.67 * ATP, 0.6]), np.array([9.9, 0.7, 0.8]), np.array([6.8])]

# dest: which state (1-based, as in MATLAB) each transition rate k{i} goes to
dest = [np.array([2, 1]), np.array([1, 1, 3]), np.array([2])]

# stepsize: whether there is a step and what the step size is for the rate k{i}
stepsize = [np.array([0, -1]), np.array([1, -1, -5]), np.array([5])]

# seqdep: whether that transition k{i} is sequence dependent
seqdep = [np.array([0, -1]), np.array([1, -1, 0]), np.array([0])]

Nsim = 100  # Nsim is the number of simulated steps

f = 0.36
r1 = 0
U = 1.9

fU = [f, U]
loop = 'TTTT'
shift_Popen = 0

Type_Output = 2  # 2, 4 column; 1, 2 column
step_limit = 0
dwell_limit = 0

# %% sequence
if Seq == 1:
    DNAhairpinseq = 'GGC TGA TAG CTG AGC GGT CGG TAT TTC AAA AGT CAA CGT ACT GAT CAC GCT GGA TCC TAG AGT CAA CGT ACT GAT CAC GCT GGA TCC TA'
    SeqPositionlimit = 15  # Popen 0.1 -> 13
elif Seq == 2:
    DNAhairpinseq = 'AAT CAG CGA TCA GAT AAC TAA CGC CCT GGG GAC TGG TAC GTC AGC TGT ATC AAG CTT CGA GAC TGG TAC GTC AGC TGT ATC AAG CTT CG'
    SeqPositionlimit = 25  # Popen 0.1 -> 21
elif Seq == 3:
    DNAhairpinseq = 'AAT CAT CGA TAA TAT AAC TAA TGC ATC GAA AAT CAG TGA AAA TCA GCT ACA ACG CCC TGG GGA CTG GTA CGT CAG CTG TAT CAA GCT TA'
    SeqPositionlimit = 35  # Popen 0.1 -> 31
elif Seq == 5:
    DNAhairpinseq = 'TTA ATA AAT AAA TAA ATA AAA TAA ATA AAG GGC GGC GGG CGG GCG GGC GGG CGG GCG GAT TAA TAA ATA AAT AAA TAA AAT AAA TAA AA'
    if F_feedback <= 9:
        SeqPositionlimit = 29
    elif F_feedback > 9:
        SeqPositionlimit = 89

# %% Popen
from ModelHairpinPopen import ModelHairpinPopen  # noqa: E402  (Zhi's code)
from SimKineticCycleSequenceTest import SimKineticCycleSequenceTest  # noqa: E402  (Zhi's code)

P_open, Pn, Hairpin_Position, Gtot = ModelHairpinPopen(1)

xpos, cdwell = SimKineticCycleSequenceTest(k, dest, stepsize, seqdep, P_open, Nsim, fU)

# %%
plt.figure()
xb, yb = stairs_xy(cdwell, xpos)
plt.step(cdwell, xpos, where='post')
x = np.asarray(cdwell, dtype=float).ravel()
y = np.asarray(xpos, dtype=float).ravel()
xs = np.repeat(x, 2)
ys = np.repeat(y, 2)
xs = xs[1:]
ys = ys[:-1]
plt.plot(xs, ys, 'r')
plt.plot(xb, yb, 'bx', label='org')

# %%
t = np.arange(0, np.max(xb) + 1e-12, 0.0001)
step_t_0 = np.concatenate([[0], xb])
y = []
for i in range(len(step_t_0) - 1):
    tx = t[(t >= step_t_0[i]) & (t <= step_t_0[i + 1])]
    y.append(np.ones(len(tx)) * yb[i])
y = np.concatenate(y)

# %%
plt.plot(t, y, 'rx', label='even data ', markersize=5)

# %%
snr = 10
yn = awgn(y, snr)
plt.plot(t, yn, 'k.', label='noisy data ')

mx = np.column_stack([t, yn])
dx = downsample(mx, 10)

plt.plot(dx[:, 0], dx[:, 1], 'r.', label='noisy data downsample')

# %%
snr = 15
yn = awgn(y, snr)

mx = np.column_stack([t, yn])
dx = downsample(mx, 10)

plt.plot(dx[:, 0], dx[:, 1], 'r.', label='less noisy data downsample')

# %%
outdir = 'D:\\Box_D\\Box\\UIUC_box_research\\progress_analysis_folder\\Recent\\meeting240501_zhisim_bnp_step\\'

np.savetxt(outdir + 'test_less_data_10.txt', dx, delimiter=',', fmt='%.5g')  # csvwrite

plt.show()
