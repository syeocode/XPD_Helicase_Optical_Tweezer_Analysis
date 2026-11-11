import numpy as np


def nowindow(num):
    win = np.ones(int(num)) / num
    return win
