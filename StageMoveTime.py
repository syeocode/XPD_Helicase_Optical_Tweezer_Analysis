import os
import numpy as np


def StageMoveTime(StartPath, stgTimeFilename):
    """
    Read _stagetime.txt files generated during time and initialize start of
    save time to 0 and outputs an array.
    _stagetime.txt saves start of data save, starts/ends moves and end of
    data save.
    Created 02/10/19 SY

    Returns (stgTime, stgPos): numpy array of times (s) and list of the
    non-time tokens (stage positions) in the file.
    """
    stgTimeFilePath = os.path.join(StartPath, stgTimeFilename)

    with open(stgTimeFilePath, 'r') as fh:
        timeArr = fh.read().split()  # textread(..., '%s')

    stgTime = []
    stgPos = []
    for token in timeArr:
        timestrarr = token.split(':')
        if len(timestrarr) > 1:
            timeHr = 60 * 60 * float(timestrarr[0])
            timeMin = 60 * float(timestrarr[1])
            timeSec = float(timestrarr[2])
            currtime = timeHr + timeMin + timeSec
            stgTime.append(currtime)
        else:
            stgPos.append(timestrarr[0])
    stgTime = np.array(stgTime, dtype=float)
    stgTime = stgTime - stgTime[0]

    # negative time resulted from over 24 hours.
    negtime_idx = stgTime < 0  # get neg time index
    stgTime[negtime_idx] = stgTime[negtime_idx] + 24 * 60 * 60

    return stgTime, stgPos
