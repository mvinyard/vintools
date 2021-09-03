
# package imports #
# --------------- #
import glob
import os

# local imports #
# ------------- #
from ...._utilities._system_utils._flexible_mkdir import _flexible_mkdir
from ._view import _view

def _get_10x_sample_name(path):

    """"""

    return path.split("/")[-3]


def _get_10x_bams(path):

    """

    Notes:
    ------

    (1) Assumes the following structure to the bam file:
        /user/specified/path/[SAMPLE]/[OUTS]/possorted_bam.bam
    """
    
    bamfile_paths = glob.glob(path + "*/outs/*.bam")
    
    SampleDict = {}

    for filepath in bamfile_paths:

        SampleDict[_get_10x_sample_name(filepath)] = filepath

    return SampleDict


def _add_outpath_to_10x_SampleDict(SampleDict, outdir):

    """Append outpath as a value associated with each key of the 10x SampleDict"""

    _flexible_mkdir(outdir)

    for [id_key, values] in SampleDict.items():
        sample_outdir = os.path.join(outdir, id_key)
        SampleDict[id_key] = [SampleDict[id_key], sample_outdir]

    return SampleDict

def _add_cleaned_outpath_to_10x_SampleDict(SampleDict):

    """Append outpath as a value associated with each key of the 10x SampleDict"""
    
    for [id_key, values] in SampleDict.items():
        cleaned_fastqs_dir = os.path.join(SampleDict[id_key][1], "cleaned_fastqs")
        cleaned_bam_dir = os.path.join(SampleDict[id_key][1], "cleaned_bam")
        SampleDict[id_key].append(cleaned_fastqs_dir)
        SampleDict[id_key].append(cleaned_bam_dir)

def _organize_10x_samples(path, outdir):

    """"""

    SampleDict = _add_outpath_to_10x_SampleDict(_get_10x_bams(path), outdir)

    return SampleDict


def _find_orphaned_samples(SampleDict, view=True):

    """"""

    OrphanedSampleDict = {}

    for [id_key, values] in SampleDict.items():

        try:
            os.path.exists(glob.glob(values[1] + "/{}*merge_region.bed".format(id_key))[0])
        except:
            OrphanedSampleDict[id_key] = values
    
    if view:
        _view(OrphanedSampleDict)
    
    return OrphanedSampleDict
