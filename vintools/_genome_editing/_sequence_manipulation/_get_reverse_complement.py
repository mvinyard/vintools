
def _get_reverse_complement(seq):

    """
    Get the reverse compliment of a DNA sequence.
    
    Parameters:
    -----------
    seq
    
    Returns:
    --------
    reverse_complement_seq
    
    Notes:
    ------
    (1) No dependencies required. Pure python. 
    """

    complement_seq = ""
    for i in seq:
        if i == "C":
            complement_seq += "G"
        elif i == "G":
            complement_seq += "C"
        elif i == "A":
            complement_seq += "T"
        elif i == "T":
            complement_seq += "A"
        elif i == "N":
            complement_seq += "N"
    reverse_complement_seq = complement_seq[::-1]
    
    return reverse_complement_seq