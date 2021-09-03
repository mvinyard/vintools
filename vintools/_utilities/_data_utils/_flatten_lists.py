def _flatten_list_of_lists(_2d_list):
    
    """Input a list of lists and output one single list"""
    
    flat_list = []
    # Iterate through the outer list
    
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list
