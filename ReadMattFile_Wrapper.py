from ReadJointFile import ReadJointFile


def ReadMattFile_Wrapper(startpath, rootfile):
    """Call this function first, it will then call the appropriate data reading
    function"""
    # before adding version numbers
    if float(rootfile[0:6]) < 90507:
        # ReadMattFile_TimeSharing was not part of the converted folder
        from ReadMattFile_TimeSharing import ReadMattFile_TimeSharing
        data = ReadMattFile_TimeSharing(startpath, rootfile)
    else:
        data = ReadJointFile(startpath, rootfile)
    return data
