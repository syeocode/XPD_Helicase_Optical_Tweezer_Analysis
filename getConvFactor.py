from matlab_helpers import Struct


def getConvFactor(StrDate):
    numDate = float(StrDate)
    ConvFactor = Struct()

    if numDate < 240410:
        ConvFactor.flconv = 123
    elif numDate > 240410:
        ConvFactor.flconv = 123  # new conv factor
    # (as in the original code, numDate == 240410 leaves flconv unset)

    if numDate < 180629:
        ConvFactor.otconvX = 891.5
        ConvFactor.otconvY = 505.7
    elif numDate < 230810:
        ConvFactor.otconvX = 883.59
        ConvFactor.otconvY = 508.93
    elif numDate > 230810:
        ConvFactor.otconvX = 896.6691
        ConvFactor.otconvY = 534.2339

    return ConvFactor
