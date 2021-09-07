
from ...._utilities._ux_utils._pystrings import _format_string_printing_font

def _view(SampleDict):
    
    print("{}: {}".format(_format_string_printing_font("Number of samples", ["BOLD"]),
                          _format_string_printing_font(str(len(SampleDict)), ["BOLD", "BLUE"])))
    print("---------------------\n")

    for [sample, values] in SampleDict.items():

        print(
            "{}: {}\n".format(
                _format_string_printing_font("Sample", ["BOLD"]),
                _format_string_printing_font(sample, ["BOLD", "RED"])))
        
        for n, val in enumerate(values):
            print("\t{}".format(_format_string_printing_font(val, ["BOLD"])))
            if n == len(values) -1:
                print("")
