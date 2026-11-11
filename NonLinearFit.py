import numpy as np
from matlab_helpers import nlinfit, nlparci


def NonLinearFit(f, spec, specFunc, init):
    """
    [fc, D, resid, sfc, sD] = NonLinearFit(f, spec, specFunc, init) returns the
    corner frequency, fc, and diffusion constant, D, from a nonlinear least
    squares fit to the function Pexp/Ptheory - 1.  This is the desired least
    squares function for spectra that have been averaged.
    Ptheory is defined by specFunc and must be defined as specFunc(f, [fc,
    D]).  init are the inital guesses for fc and D, [fc, D].

    Jeff Moffitt
    8/5/06
    """
    spec = np.asarray(spec, dtype=float)
    f = np.asarray(f, dtype=float)

    # Define Pexp/Ptheory for the fitting routine
    func = lambda para, f: spec / specFunc(f, para)

    # Call the nonlinear least squares minimization routine
    paraFit, resid, J = nlinfit(f, np.ones(len(spec)), func, init)

    # Assign fit values
    fc = paraFit[0]
    D = paraFit[1]

    # Calculate 95% Confidence Intervals for the fit
    ci = nlparci(np.real(paraFit), np.real(resid), np.real(J))  # Kludge!!

    # Calculate 1 sigma from these intervals
    sfc = (ci[0, 1] - ci[0, 0]) / 4
    sD = (ci[1, 1] - ci[1, 0]) / 4

    return fc, D, resid, sfc, sD
