import glob, os
import numpy as np

def _format_search_within_path(path):
    if path.endswith("*"):
        return path
    elif path.endswith("/"):
        return path + "*"
    else:
        return path + "/*"

def _identify_sample_path_structure(path, file_extension):
    
    if path.count("*") > 1:
        path = path.strip("*" + file_extension)
    
    stripped_path = np.array(path.split("/"))
    wildcard_loc = np.where(stripped_path == '*')[0] + 1

    return wildcard_loc.item()

def _breakdown_single_path(sample_path, sample_dir_level, file_extension):

    path_components = []
    for path_component in sample_path.split("/"):
        path_components.append(path_component)
    
    sample_name = os.path.basename("/".join(path_components[:sample_dir_level]))
        
    return sample_name

def _find_sample_files(path, file_extension, sample_dir_level=False, silent=False):
    
    """Given a directory and file extension, finds files of the extension type and matches them to samples."""
    
    FileDict = {}
    sample_paths = glob.glob(_format_search_within_path(path) + file_extension)
    
    if not sample_dir_level:
        sample_dir_level = _identify_sample_path_structure(path, file_extension)

    if not silent:
        print("There were {} {} paths detected.\n".format(str(len(sample_paths)), "." + file_extension))
        print("Given the current configuration, the sample names would be:\n")
        
    for path in sample_paths:
        sample_name = _breakdown_single_path(path, sample_dir_level, file_extension)
        if not silent:
            print("\t", sample_name)
        FileDict[sample_name] = path
    
    return FileDict

def _plotCorrelation(title, results_file="multiBigwigsummary.results.npz", silent=False):
    
    """Takes as input, the results.npz file generated in the previous step."""
    
    for correlation_method in ['spearman', 'pearson']:
        for plot_method in ["heatmap", "scatterplot"]:
            plotfile = ".".join([title, correlation_method, plot_method, "correlation.png"])
            correlation_matrix_file = ".".join([title, correlation_method, plot_method, "correlation_matrix.tsv"])
            os.system("plotCorrelation --corData {} --corMethod {} --whatToPlot {} --plotFile {} --outFileCorMatrix {}".format(results_file, correlation_method, plot_method, plotfile, correlation_matrix_file))
            if not silent:
                print(plotfile)