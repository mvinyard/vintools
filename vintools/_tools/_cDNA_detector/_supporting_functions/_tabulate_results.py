from ...._plotting._annotate_legend import _annotate_legend
from ...._plotting.color_palettes import vin_colors
from ...._plotting._modify_ax_spines import _modify_ax_spines
from ...._utilities._ux_utils._pystrings import _format_string_printing_font

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import glob
import os


def _fetch_gene_statistics_source_filtered(results_dir, sample_id):

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
    #     print(results_dir + "{}/*gene_statistics.filtered.source_filtered.tsv".format(sample_id))

    filtered_gene_stats_path = glob.glob(
        results_dir
        + "{}/*gene_statistics.filtered.source_filtered.tsv".format(sample_id)
    )[0]
    #     print(filtered_gene_stats_path)
    filtered_gene_stats = pd.read_csv(filtered_gene_stats_path, sep="\t")

    return filtered_gene_stats


def _fetch_merged_regions(results_dir, sample_id):

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

    #     print(results_dir + "{}/*merge_region.tsv".format(sample_id))
    merged_region_path = glob.glob(
        results_dir + "{}/*merge_region.tsv".format(sample_id)
    )[0]
    #     print(merged_region_path)
    merged_regions = pd.read_csv(merged_region_path, sep="\t")
    return merged_regions


def _echo_SampleDict_composition(ResultsDict):

    print("\nSamples:")
    print("--------\n")

    for key in ResultsDict.keys():
        try:
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

        except:
            pass


def _fetch_results(results_dir, silent=False):

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
    ResultsDict = {}
    missing_results = []

    for n, sample in enumerate(samples):
        sample_id = os.path.basename(sample)
        try:
            ResultsDict[sample_id] = {}

            ResultsDict[sample_id][
                "gene_statistics_source_filtered"
            ] = _fetch_gene_statistics_source_filtered(results_dir, sample_id)

            ResultsDict[sample_id]["merged_regions"] = _fetch_merged_regions(
                results_dir, sample_id
            )
        except:
            missing_results.append(sample_id)
            del ResultsDict[sample_id]
    if not silent:
        _echo_SampleDict_composition(ResultsDict)

    return ResultsDict, missing_results


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
        try:
            n_fragments_passing.append(tabulated_results[key]["n_fragments_passing"])
            n_fragments.append(tabulated_results[key]["n_fragments"])
        except:
            pass
    return n_fragments, n_fragments_passing


def _plot_tabulated_results(tabulated_results):

    cols = vin_colors()

    fig = plt.figure(figsize=(14, 12))
    gridspec = GridSpec(1, 1)
    bar_ax = fig.add_subplot(gridspec[0, 0])
    bar_ax_spines = _modify_ax_spines(bar_ax)
    bar_ax_spines.delete(["top", "right", "left"])

    n_fragments, n_fragments_passing = _grab_fragments(tabulated_results)

    bar_ax.bar(
        np.array(range(len(n_fragments))) - 0.25,
        n_fragments,
        width=0.4,
        color=cols[3],
        label="Num fragments",
    )
    bar_ax.bar(
        np.array(range(len(n_fragments_passing))) + 0.25,
        n_fragments_passing,
        width=0.4,
        color=cols[9],
        label="Num passing fragments",
    )
    x = bar_ax.set_xticks(np.array(range(0, len([*tabulated_results.keys()]))))
    x = bar_ax.set_xticklabels([*tabulated_results.keys()], rotation=40, ha="right")

    max_y = np.max([np.max(n_fragments_passing), np.max(n_fragments)])
    y_ticks = np.round(np.linspace(0, max_y, 8))
    y = bar_ax.set_yticks(y_ticks)
    _annotate_legend(bar_ax, loc=1)


def _tabulate_results(results_dir, plot=True, silent=False):

    ResultsDict, missing_results = _fetch_results(results_dir, silent=silent)

    tabulated_results = {}

    for key in ResultsDict.keys():
        try:

            tabulated_results[key] = {}
            merged_regions_df = ResultsDict[key]["merged_regions"]
            tabulated_results[key]["n_fragments"] = merged_regions_df.shape[0]
            tabulated_results[key]["n_fragments_passing"] = merged_regions_df.loc[
                merged_regions_df["filter"] == "pass"
            ].shape[0]
        except:
            pass

    if plot:
        _plot_tabulated_results(tabulated_results)

    return tabulated_results, missing_results
