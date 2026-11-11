from matlab_helpers import Struct, downsample, fieldnames


def downsample_trapdata(trapdatai, factor):
    """
    downsample all data after trap_data read
    190310 SY
    (All fields except the last one are downsampled; the last field is
    assumed to be the sampling period and is multiplied by factor.)
    """
    fieldn = fieldnames(trapdatai)
    trapdata = Struct()
    for name in fieldn[:-1]:
        data = trapdatai[name]
        data = downsample(data, factor)
        trapdata[name] = data
    name = fieldn[-1]
    sampp = trapdatai[name]
    trapdata[name] = sampp * factor
    return trapdata
