import numpy as np
from matlab_helpers import inputdlg


def getnumbers(title, prompts, defaults):
    """Ask for numbers (console version of the MATLAB input dialog)."""
    cells = inputdlg(prompts, title, defaults)
    N = max(np.shape(prompts)) if len(prompts) else 0
    numbers = np.zeros(N)
    for i in range(N):
        numbers[i] = float(eval(cells[i], {"__builtins__": {}}, {}))  # like str2num
    return numbers
