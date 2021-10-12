def _read_txt(filepath):

    """
    Read a text file given a filepath.
    
    Parameters:
    -----------
    filepath
        path to textfile
        type:str
    
    Returns:
    --------
    lines
        lines of content in the file appended into a list
        type: list
    
    Notes:
    ------
    (1) Skips lines without content. 
    """

    lines = []

    f = open(filepath, "r")
    for line in f.readlines():
        if len(line) == 1:
            continue
        else:
            lines.append(line.strip("\n"))

    return lines