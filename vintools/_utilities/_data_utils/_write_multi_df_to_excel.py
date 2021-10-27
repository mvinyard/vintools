
# _write_multi_df_to_excel.py
__module_name__ = "_write_multi_df_to_excel.py"
__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# package imports #
# --------------- #
import pandas as pd

# local imports #
# ------------- #
from ._check_fix_file_extension import _check_fix_file_extension


def _write_multi_df_to_excel(
    df_list, workbook_path="./new_workbook.xlsx", sheetnames=False, index=False, silent=False
):
    
    """
    Write one or more pandas.DataFrames to individual sheets of an Excel Workbook (workbook.xlsx). 
    
    Parameters:
    -----------
    df_list
        type: list(pd.DataFrame, ...)
    
    workbook_path
        type: str
    
    sheetnames
    If passed, Must be equal in length to df_list
        default: False
        type: bool, or if passed: list(str, ...)
    index
        default: False
        type: bool
    silent
        default: False
        type: bool
        
    Returns:
    --------
    None
        writes to excel workbook (.xlsx)
        
    Notes:
    ------
    (1) if 
    """

    workbook_path = _check_fix_file_extension(workbook_path, ".xlsx", silent)

    if sheetnames:
        assert len(df_list) == len(sheetnames), print(
            "Length of df_list must equal the length of sheet names passed."
        )
    with pd.ExcelWriter(workbook_path) as writer:

        if not silent:
            print("Writing to: {}\n".format(workbook_path))

        for n, df_ in enumerate(df_list):
            if not sheetnames:
                sheet = "Sheet_{}".format(n)
            else:
                sheet = sheetnames[n]
            if not silent:
                print("\tSheet {}:\t{}".format(str(int(n + 1)), sheet))
            df_.to_excel(writer, sheet, index=index)