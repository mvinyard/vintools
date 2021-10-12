def _handle_kwargs(verbose=True, **kwargs):

    """
    Take input keyword arguments passed into a wrapping function and organize into a callable dictionary.

    Parameters:
    -----------
    **kwargs
        keyword arguments followed by assigned values of any type

    Returns:
    --------
    ArgDict
        type: dict

    Notes:
    ------
    """

    KwargDict = {}

    for key, value in kwargs.items():
        if verbose:
            _print_input_kwarg(key, value)
        KwargDict[key] = value

    return KwargDict
