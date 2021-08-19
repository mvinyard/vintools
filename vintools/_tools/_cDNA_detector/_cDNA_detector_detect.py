
# package imports #
# --------------- #
import glob
import os

# local imports #
# ------------- #
from ..._utilities._flexible_mkdir import _flexible_mkdir

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


def _format_outdir(self, outdir):

    if not outdir:
        self.outdir = "./"
        print(
            "Outdir not defined. Will return results to current working dir: {}".format(
                os.getcwd()
            )
        )
    else:
        self.outdir = outdir


def _detect(
    script, bam, sample_id, gene_model, n_cores, cDNA_detector_outdir_sample, dry_run
):

    """"""

    print("Now runnning cDNA Detector detect module on sample: {}".format(sample_id))

    cDNA_detector_detect_command = "{} detect --bam {} --sample_id {} --gene_model {} --n_threads {} --output_dir {}\n".format(
        script, bam, sample_id, gene_model, n_cores, cDNA_detector_outdir_sample
    )
    if dry_run:
        print(cDNA_detector_detect_command)
    else:
        os.system(cDNA_detector_detect_command)