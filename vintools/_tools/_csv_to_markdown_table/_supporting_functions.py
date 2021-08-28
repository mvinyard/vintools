import numpy as np

def _add_bookends(row):
    return "".join(["|", "|".join(row), "|"])


def _make_header(df):
    return _add_bookends(df.columns) + "\n"


def _make_header_delimiter(ncols):
    return "".join(np.append(np.repeat(["|--"], ncols), np.array(["|"]))) + "\n"


def _make_TableDict(df):

    table = {}
    table["header"] = _make_header(df)
    table["header_delimiter"] = _make_header_delimiter(len(df.columns))

    for idx, row in df.iterrows():
        table[idx] = _add_bookends(row.values.astype(str)) + "\n"
    return table

def _view_formatted_table(table):
    
    for [key, row] in table.items():
        print(row.strip("\n"))
        
def _write_dict_to_markdown_table(TableDict, outfile):

    """"""
    f = open(outfile, "w")
    for [key, val] in TableDict.items():
        f.write(val)
    f.close()