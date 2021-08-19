
# package imports #
# --------------- #
import glob
import os

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

    bam_path = path + "*/outs/*.bam"

    bamfile_paths = glob.glob(bam_path)

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


def _organize_10x_samples(path, outdir):

    """"""

    SampleDict = _add_outpath_to_10x_SampleDict(_get_10x_bams(path), outdir)

    return SampleDict