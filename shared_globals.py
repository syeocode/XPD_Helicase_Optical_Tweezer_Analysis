"""
shared_globals.py

Python replacement for the MATLAB `global` variables used by the analysis
code (WLC_param, construct, datadirectories, calparams, ...) and for the
MATLAB "base workspace" used by assignin/evalin.

Set them once at the start of a session, e.g.

    import shared_globals as G
    from matlab_helpers import Struct
    G.WLC_param = Struct(hds=0.34, hss=0.59, Pds=50, Pss=1.0,
                         Sds=1200, Sss=800, kT=4.14)
    G.construct = '5Fork-Hairpin-10dT'
    G.datadirectories = 'E:\\data\\'

and every converted function that used `global X` in MATLAB reads `G.X`.
"""

# Worm-like-chain parameters (struct with hds, hss, Pds, Pss, Sds, Sss, kT)
WLC_param = None
# Name of the DNA construct (see ConstructConstants)
construct = None

# Data folders
datadirectories = ''
datadirectory = ''
XPDdatadir = ''
savedirectory = ''
analysis = ''
fbslash = '\\'

# Calibration
calparams = None          # [beadA beadB fsamp fXYhigh fSumHigh fLow avwin]
mostrecentcal = 0
trap_computer = 0

# Stand-in for MATLAB's base workspace (assignin / evalin('base', ...))
base_workspace = {}
