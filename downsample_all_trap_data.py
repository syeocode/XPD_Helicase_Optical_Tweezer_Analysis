import copy
from matlab_helpers import downsample


def downsample_all_trap_data(data, factor):
    """
    downsample all data associated with the trap and update sampling period
    value.  Most appropriately called right after opening the data.  Leaves
    fluorescence data alone.  Also only changes for 2 traps (A, B, ignores
    possiblility of C).
    111202 mjc
    """
    data = copy.copy(data)  # MATLAB passes structs by value
    data.sampperiod = data.sampperiod * factor

    data.A_Y = downsample(data.A_Y, factor)
    data.A_X = downsample(data.A_X, factor)
    data.B_Y = downsample(data.B_Y, factor)
    data.B_X = downsample(data.B_X, factor)

    data.A_Sum = downsample(data.A_Sum, factor)
    data.B_Sum = downsample(data.B_Sum, factor)

    data.trappos1 = downsample(data.trappos1, factor)
    data.trappos2 = downsample(data.trappos2, factor)

    return data
