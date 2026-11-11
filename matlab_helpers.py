"""
matlab_helpers.py

Small Python stand-ins for MATLAB built-ins that the converted analysis code
relies on (struct, std, smooth, movmean, downsample, interp1, nlinfit/nlparci,
fminsearch, findpeaks, inputdlg, ...). Each one mimics MATLAB's behaviour so the
converted functions give the same numbers as the original .m files.

Created during the MATLAB -> Python conversion of New_Default_Analysis.
"""
import numpy as np
from types import SimpleNamespace


# ---------------------------------------------------------------------------
# struct
# ---------------------------------------------------------------------------
class Struct(SimpleNamespace):
    """MATLAB-like struct: access fields with s.name, s['name'] or getattr.

    Fields keep the order in which they were created (like MATLAB's
    fieldnames). Use isfield(s, 'name') to test for a field.
    """

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __contains__(self, key):
        return hasattr(self, key)

    def fieldnames(self):
        return list(vars(self).keys())

    def copy(self):
        return Struct(**vars(self))


def isfield(s, name):
    if isinstance(s, dict):
        return name in s
    return hasattr(s, name)


def fieldnames(s):
    if isinstance(s, dict):
        return list(s.keys())
    return list(vars(s).keys())


def struct_from_dict(d):
    """Recursively turn a dict (e.g. from scipy.io.loadmat) into a Struct."""
    s = Struct()
    for k, v in d.items():
        if k.startswith('__'):
            continue
        s[k] = struct_from_dict(v) if isinstance(v, dict) else v
    return s


# ---------------------------------------------------------------------------
# statistics / array helpers
# ---------------------------------------------------------------------------
def mstd(x, axis=None):
    """MATLAB std (normalised by N-1; std of a single value is 0)."""
    x = np.asarray(x, dtype=float)
    if axis is None and x.size == 1:
        return 0.0
    return np.std(x, axis=axis, ddof=1)


def isempty(x):
    if x is None:
        return True
    if isinstance(x, (str, list, tuple, dict)):
        return len(x) == 0
    return np.size(x) == 0


def downsample(x, n):
    """MATLAB downsample: keep every n-th sample starting with the first.
    For 2-D input, rows are downsampled (MATLAB operates on columns)."""
    x = np.asarray(x)
    return x[::int(n)]


def movmean(x, k):
    """MATLAB movmean(x, k) for a vector, with shrinking windows at the ends.
    Odd k: centred window. Even k: k/2 points before and k/2-1 after."""
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    k = int(k)
    if k % 2 == 1:
        before = after = (k - 1) // 2
    else:
        before = k // 2
        after = k // 2 - 1
    c = np.concatenate(([0.0], np.cumsum(x)))
    idx = np.arange(n)
    lo = np.clip(idx - before, 0, n)
    hi = np.clip(idx + after + 1, 0, n)
    return (c[hi] - c[lo]) / (hi - lo)


def movmedian(x, k):
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    k = int(k)
    if k % 2 == 1:
        before = after = (k - 1) // 2
    else:
        before = k // 2
        after = k // 2 - 1
    return np.array([np.median(x[max(0, i - before):min(n, i + after + 1)]) for i in range(n)])


def smooth(y, span=5, method='moving'):
    """MATLAB smooth(y, span, 'moving'): moving average whose window shrinks
    near the ends so that it stays symmetric (span forced to be odd)."""
    if method != 'moving':
        raise NotImplementedError("Only the 'moving' method is implemented")
    y = np.asarray(y, dtype=float).ravel()
    n = y.size
    span = int(min(span, n))
    if span % 2 == 0:
        span -= 1
    if span <= 1:
        return y.copy()
    half = (span - 1) // 2
    c = np.concatenate(([0.0], np.cumsum(y)))
    idx = np.arange(n)
    h = np.minimum(np.minimum(idx, n - 1 - idx), half)
    return (c[idx + h + 1] - c[idx - h]) / (2 * h + 1)


def interp1(x, y, xi):
    """MATLAB interp1 (linear); points outside the data range give NaN."""
    from scipy.interpolate import interp1d
    f = interp1d(np.asarray(x, float), np.asarray(y, float), kind='linear',
                 bounds_error=False, fill_value=np.nan, assume_sorted=False)
    return f(np.asarray(xi, float))


def lfilter_ma(N, x):
    """MATLAB filter(ones(1,N), N, x)."""
    from scipy.signal import lfilter
    return lfilter(np.ones(int(N)), float(N), np.asarray(x, float))


def findpeaks(y):
    """MATLAB [pks, locs] = findpeaks(y) (0-based locs)."""
    from scipy.signal import find_peaks
    y = np.asarray(y, float).ravel()
    locs, _ = find_peaks(y)
    return y[locs], locs


def logical_index(x, mask):
    """MATLAB x(mask) where mask may be shorter than x (missing = false)."""
    x = np.asarray(x)
    mask = np.asarray(mask, dtype=bool)
    return x[:mask.size][mask]


# ---------------------------------------------------------------------------
# optimisation
# ---------------------------------------------------------------------------
def fminsearch(fun, x0):
    """MATLAB fminsearch (Nelder-Mead with MATLAB's default tolerances)."""
    from scipy.optimize import minimize
    x0 = np.atleast_1d(np.asarray(x0, dtype=float))
    res = minimize(lambda v: fun(v[0] if x0.size == 1 else v), x0,
                   method='Nelder-Mead',
                   options={'xatol': 1e-4, 'fatol': 1e-4,
                            'maxiter': 200 * x0.size, 'maxfev': 200 * x0.size})
    return res.x[0] if x0.size == 1 else res.x


def nlinfit(X, y, modelfun, beta0):
    """MATLAB [beta, r, J] = nlinfit(X, y, modelfun, beta0)
    (ordinary least squares, Levenberg-Marquardt).
    modelfun(beta, X) -> predicted y.  r = y - f,  J = d f / d beta."""
    from scipy.optimize import least_squares
    y = np.asarray(y, dtype=float)
    beta0 = np.asarray(beta0, dtype=float)

    def resid(b):
        return np.real(y - modelfun(b, X))

    res = least_squares(resid, beta0, method='lm', x_scale='jac',
                        xtol=1e-8, ftol=1e-8, max_nfev=2000)
    beta = res.x
    r = resid(beta)
    J = -res.jac  # jacobian of the model, as MATLAB returns it
    return beta, r, J


def nlparci(beta, resid, J, alpha=0.05):
    """MATLAB ci = nlparci(beta, resid, 'jacobian', J) (95% by default)."""
    from scipy.stats import t
    beta = np.asarray(beta, dtype=float)
    resid = np.asarray(resid, dtype=float).ravel()
    J = np.asarray(J, dtype=float)
    n, p = J.shape
    dof = n - p
    mse = np.sum(resid ** 2) / dof
    # same as MATLAB: use the R factor of J
    _, R = np.linalg.qr(J)
    Rinv = np.linalg.pinv(R)
    se = np.sqrt(np.sum(Rinv ** 2, axis=1) * mse)
    delta = se * t.ppf(1 - alpha / 2, dof)
    return np.column_stack([beta - delta, beta + delta])


# ---------------------------------------------------------------------------
# user input (replacement for inputdlg)
# ---------------------------------------------------------------------------
def inputdlg(prompts, title='', defaults=None):
    """Console replacement for MATLAB inputdlg: returns a list of strings.
    Press Enter to accept the default shown in brackets."""
    if isinstance(prompts, str):
        prompts = [prompts]
    if defaults is None:
        defaults = [''] * len(prompts)
    if title:
        print(title)
    out = []
    for p, d in zip(prompts, defaults):
        ans = input(f'{p} [{d}]: ')
        out.append(ans if ans.strip() != '' else str(d))
    return out


# ---------------------------------------------------------------------------
# file helpers
# ---------------------------------------------------------------------------
def fread_f64(fid, n=1):
    """fread(fid, n, 'float64') on a big-endian file."""
    v = np.fromfile(fid, dtype='>f8', count=n)
    return float(v[0]) if n == 1 else v


def loadmat_struct(path):
    """load('file.mat') -> Struct of the variables in the file."""
    from scipy.io import loadmat
    d = loadmat(path, squeeze_me=True, struct_as_record=False)
    return struct_from_dict(d)
