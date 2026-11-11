from matlab_helpers import inputdlg


def getnumbers_fext(title, prompts, defaults):
    """Ask for force-extension parameters; returns the 8 answers as strings."""
    cells = inputdlg(prompts, title, defaults)

    fext_array = [None] * 8
    fext_array[0] = cells[0]
    fext_array[1] = cells[1]  # [nm/bp] interphosphate distance
    fext_array[2] = cells[2]  # nm/nt
    fext_array[3] = cells[3]  # nm
    fext_array[4] = cells[4]  # nm
    fext_array[5] = cells[5]  # pN
    fext_array[6] = cells[6]  # pN
    fext_array[7] = cells[7]  # Zhi's values (dsDNA)
    return fext_array
