
from ...._plotting._annotate_legend import _annotate_legend
from ...._plotting.color_palettes import vin_colors
from ...._plotting._modify_ax_spines import _modify_ax_spines
from ...._utilities._pystrings import _format_string_printing_font

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os

def _fetch_gene_statistics_source_filtered(results_dir):

    """
    Gets: sample.gene_statistics.filtered.source_filtered.tsv

    Parameters:
    -----------
    SampleDict

    results_dir

    Returns:
    --------
    SampleDict

    Notes:
    ------

    """

    filtered_gene_stats_paths = glob.glob(
        results_dir + "*/*gene_statistics.filtered.source_filtered.tsv"
    )

    filtered_gene_stats = []
    for path in filtered_gene_stats_paths:
        filtered_gene_stats.append(pd.read_csv(path, sep="\t"))

    return filtered_gene_stats


def _fetch_merged_regions(results_dir):

    """
    Gets: sample.gene_statistics.filtered.source_filtered.tsv

    Parameters:
    -----------
    results_dir

    Returns:
    --------

    Notes:
    ------

    """

    merged_region_paths = glob.glob(results_dir + "*/*merge_region.tsv")

    merged_regions = []
    for path in merged_region_paths:
        merged_regions.append(pd.read_csv(path, sep="\t"))

    return merged_regions


def _echo_SampleDict_composition(ResultsDict):

    print("\nSamples:")
    print("--------\n")

    for key in ResultsDict.keys():
        print("{}".format(_format_string_printing_font(key, "BOLD")))
        keylist = [*ResultsDict[key].keys()]
        print("\tNum genes: {}".format(_print_n_genes(ResultsDict, key), "BOLD"))
        print(
            "\tNum passing fragments: {}".format(
                _print_n_genes(ResultsDict, key), "BOLD"
            )
        )
        print("\n\t{}".format(_format_string_printing_font("DataFrames:", "BOLD")))
        for subkey in keylist:
            ResultsDict[key][subkey]
            print("\t  {}".format(subkey))
            
            
def _fetch_results(results_dir):

    """
    Get gene statistics source filtered and merged region. Format as pandas dfs.

    Parameters:
    -----------
    results_dir

    Results:
    --------
    SampleDict

    Notes:
    ------
    """

    samples = glob.glob(results_dir + "*")

    gene_statistics_source_filtered = _fetch_gene_statistics_source_filtered(
        results_dir
    )
    merged_regions = _fetch_merged_regions(results_dir)

    ResultsDict = {}

    for n, sample in enumerate(samples):
        sample_id = os.path.basename(sample)
        try:
            ResultsDict[sample_id] = {}
            ResultsDict[sample_id][
                "gene_statistics_source_filtered"
            ] = gene_statistics_source_filtered[n]
            ResultsDict[sample_id]["merged_regions"] = merged_regions[n]
        except:
            ResultsDict[sample_id] = None
            
    _echo_SampleDict_composition(ResultsDict)
    
    return ResultsDict


def _print_n_genes(ResultsDict, key):
    return ResultsDict[key]["gene_statistics_source_filtered"]["gene_id"].nunique()


def _print_n_passing_fragments(ResultsDict, key):
    df = ResultsDict[key]["merged_regions"]
    n_passing_fragments = df.loc[df.filter == "pass"].shape[0]
    return n_passing_fragments

def _grab_fragments(tabulated_results):

    n_fragments = []
    n_fragments_passing = []
    for key in tabulated_results.keys():
        n_fragments_passing.append(tabulated_results[key]["n_fragments_passing"])
        n_fragments.append(tabulated_results[key]["n_fragments"])
        
    return n_fragments, n_fragments_passing

def _plot_tabulated_results(tabulated_results):

    cols = vin_colors(separable=True)

    fig = plt.figure(figsize=(7, 6))
    gridspec = GridSpec(1, 1)
    bar_ax = fig.add_subplot(gridspec[0, 0])
    bar_ax_spines = _modify_ax_spines(bar_ax)
    bar_ax_spines.delete(["top", "right", "left"])
    
    n_fragments, n_fragments_passing = _grab_fragments(tabulated_results)

    bar_ax.bar(
        np.array(range(len(n_fragments))) - 0.25,
        n_fragments,
        width=0.4,
        color=cols[0],
        label="Num fragments",
    )
    bar_ax.bar(
        np.array(range(len(n_fragments_passing))) + 0.25,
        n_fragments_passing,
        width=0.4,
        color=cols[1],
        label="Num passing fragments",
    )
    x = bar_ax.set_xticks(range(0, 5))
    x = bar_ax.set_xticklabels([*tabulated_results.keys()])
    
    max_y = np.max([np.max(n_fragments_passing), np.max(n_fragments)])
    y_ticks = np.round(np.linspace(0, max_y, 8))
    y = bar_ax.set_yticks(y_ticks)
    _annotate_legend(bar_ax, loc=1)
    
    
def _tabulate_results(results_dir):
    
    ResultsDict = _fetch_results(results_dir)

    tabulated_results = {}

    for key in ResultsDict.keys():

        tabulated_results[key] = {}
        merged_regions_df = ResultsDict[key]["merged_regions"]
        tabulated_results[key]["n_fragments"] = merged_regions_df.shape[0]
        tabulated_results[key]["n_fragments_passing"] = merged_regions_df.loc[
            merged_regions_df["filter"] == "pass"
        ].shape[0]
    
    _plot_tabulated_results(tabulated_results)
    
    return tabulated_results