import glob, os


def _create_10x_sample_OutDict(outsdir_path):

    outfilelist = glob.glob(os.path.join(outsdir_path, "outs/*"))

    OutDict = {}
    for file in outfilelist:
        OutDict[os.path.basename(file)] = file

    return OutDict


def _generate_10x_cellranger_ResultsDict(path_to_results, samples):

    paths_in_rundir = glob.glob(path_to_results)

    ResultsDict = {}
    ResultsDict["run_path"] = path_to_results
    ResultsDict["samples"] = {}

    for path in paths_in_rundir:
        sample = os.path.basename(path)
        if (os.path.isdir(path)) & (sample in samples):
            ResultsDict["samples"][sample] = {}
            ResultsDict["samples"][sample]["results_path"] = path
            ResultsDict["samples"][sample]["outs"] = _create_10x_sample_OutDict(
                outsdir_path=path
            )

    return ResultsDict
