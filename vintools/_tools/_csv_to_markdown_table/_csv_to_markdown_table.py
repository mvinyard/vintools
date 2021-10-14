# local imports #
# ------------- #
from ..._utilities._ux_utils._pystrings import _format_string_printing_font
from ._supporting_functions import (
    _make_TableDict,
    _write_dict_to_markdown_table,
    _view_formatted_table,
)

# package imports #
# --------------- #
import pandas as pd


class create_md_table:
    def __init__(self, outfile=False):

        """"""

        if outfile:
            self.outfile = outfile

    def read_df(self, infile):

        """
        Reads .csv file using pandas to self.df
        
        Parameters:
        -----------
        infile
            path to csv

        Returns:
        --------
        None
            self.df is created. 

        Notes:
        ------    
        """

        self.infile = infile

        self.df = pd.read_csv(self.infile)

    def make_table(self,):

        """
        reformats pandas DataFrame to a dictionary as an intermediate to writing a markdown table. 
        
        Parameters:
        -----------
        None

        Returns:
        --------
        None
            self.table is created. 

        Notes:
        ------
        This intermediate step is mostly so the user can view the contents of the table before saving the dict to a file. 
        """

        self.table = _make_TableDict(self.df)

    def view(self):

        """
        Preview the formatted markdown table. 
        
        Parameters:
        -----------
        None

        Returns:
        --------
        None
            prints a formatted view of self.table. 

        Notes:
        ------
        """

        _view_formatted_table(self.table)

    def write(self, outfile=False):

        """
        write reformated dictionary to markdown table file.
        
        Parameters:
        -----------
        None

        Returns:
        --------
        None
            self.table is created. 

        Notes:
        ------    
        """

        if outfile:
            self.outfile = outfile
        else:
            print("No outfile path defined!")

        _write_dict_to_markdown_table(self.table, self.outfile)

        print(
            "{} converted to markdown table and saved: {}".format(
                _format_string_printing_font(self.infile, ["BOLD", "RED"]),
                _format_string_printing_font(self.outfile, ["BOLD", "BLUE"]),
            )
        )


def quicktable(infile, outfile="table.md"):

    """
    Parameters:
    -----------
    infile
        path to csv
        
    outfile
        path to written markdown table
    
    Returns:
    --------
    None
        saves markdown file
    
    Notes:
    ------    
    """

    MarkdownTable = create_md_table()
    MarkdownTable.read_df(infile=infile)
    MarkdownTable.make_table()
    MarkdownTable.write(outfile=outfile)
