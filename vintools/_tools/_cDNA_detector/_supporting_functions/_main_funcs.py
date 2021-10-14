# package imports #
# --------------- #
import os
import glob

# local imports #
# ------------- #
from ...._utilities._ux_utils._pystrings import _format_string_printing_font
from ...._utilities._system_utils._clone_GitHub_repo import _clone_GitHub_repo
from ...._utilities._system_utils._flexible_mkdir import _flexible_mkdir

from ._make_10x_SampleDict import (
    _add_cleaned_outpath_to_10x_SampleDict,
    _find_orphaned_samples,
)


def _setup_cDNA_detector(cDNA_detector_repo_path, script):

    """"""

    print(
        "{} repository not found... cloning locally now.\n".format(
            _format_string_printing_font("cDNA-detector", ["BOLD", "BLUE"])
        )
    )
    _clone_GitHub_repo(
        repository="https://github.com/rheinbaylab/cDNA-detector.git",
        destination=cDNA_detector_repo_path,
    )

    os.system("chmod a+x {}".format(script))

    # install dependencies
    os.system(
        "pip install --user pandas numpy scipy pysam statsmodels biopython gffpandas joblib"
    )
    os.system("conda install -c bioconda blast==2.9.0 -y")

    print(
        "{} package and dependencies installed!".format(
            _format_string_printing_font("cDNA-detector", ["BOLD", "BLUE"])
        )
    )


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
    script,
    bam,
    sample_id,
    gene_model,
    n_cores,
    fast,
    cDNA_detector_outdir_sample,
    dry_run,
):

    """"""

    print(
        "Now runnning cDNA Detector {} module on sample: {}\n".format(
            _format_string_printing_font("detect", ["BOLD"]),
            _format_string_printing_font(sample_id, ["BOLD"]),
        )
    )
    if fast:
        fast = "&"
    else:
        fast = ""

    cDNA_detector_detect_command = "{} detect --bam {} --sample_id {} --gene_model {} --n_threads {} --output_dir {} {}\n".format(
        script, bam, sample_id, gene_model, n_cores, cDNA_detector_outdir_sample, fast
    )
    if dry_run:
        print(cDNA_detector_detect_command)
    else:
        os.system(cDNA_detector_detect_command)


def _prepend_gene_model(prepare_outdir):

    """"""

    gene_model = glob.glob(prepare_outdir + "/*.saf")[0]

    print(
        "\nUse the following gene model for `{}`:\n\n\t{}".format(
            _format_string_printing_font("cDNA-detector detect", ["BOLD"]),
            _format_string_printing_font(gene_model, ["BOLD", "BLUE"]),
        )
    )

    return gene_model


def _prepare(script, ref_10x):

    """"""

    annotations = os.path.join(ref_10x, "genes/genes.gtf")
    genome = os.path.join(ref_10x, "fasta/genome.fa")
    gene_model_dir = os.path.join(ref_10x, "gene_model")
    _flexible_mkdir(gene_model_dir)

    print("Preparing cDNA-detector gene model... this should take 5-7 minutes...")

    executable = "{} prepare --annotation {} --genome {} --output_dir {}".format(
        script, annotations, genome, gene_model_dir
    )
    os.system(executable)

    assert os.path.exists(gene_model_dir), print(
        "The gene model was not properly created, please investigate and re-run."
    )

    return _prepend_gene_model(gene_model_dir)


def _download_bam2fastq_10x_binary(
    binary_https="https://github.com/10XGenomics/bamtofastq/releases/download/v1.3.5/bamtofastq_linux",
):

    """"""

    print(
        "Path to binary for bam2fastq_linux was not found... downloading. No installation required."
    )

    binary_path_destination = os.path.dirname(_get_pypi_package_loc())
    binary = os.path.join(binary_path_destination, os.path.basename(binary_https))

    os.system("wget {} -P {}".format(binary_https, binary_path_destination))
    os.system("chmod 700 {}".format(binary))
    return binary


def _bam2fastq_10x(bam2fastq_binary, bam_file, fastq_outpath, n_threads):

    """"""

    if not os.path.exists(bam2fastq_binary):
        bam2fastq_path = _download_bam2fastq_10x_binary()

    executable = "{} --nthreads {} {} {}".format(
        bam2fastq_binary, n_threads, bam_file, fastq_outpath
    )
    os.system(executable)


def _run_clean_preflight(self):

    """"""

    _add_cleaned_outpath_to_10x_SampleDict(self.SampleDict)
    try:
        self.OrphanedSampleDict
    except:
        self.OrphanedSampleDict = _find_orphaned_samples(self.SampleDict, view=False)

    self.clean_preflight_complete = True


def _clean(script, bam, sample_id, region, output_dir):

    """"""

    _flexible_mkdir(output_dir)
    executable = "{} clean --bam {} --sample_id {} --region {} --output_dir {}".format(
        script, bam, sample_id, region, output_dir
    )

    os.system(executable)


def _clean_bam_batch(script, SampleDict, OrphanedSampleDict):

    for [id_key, values] in SampleDict.items():
        if not id_key in OrphanedSampleDict.keys():

            print(
                "Now runnning cDNA Detector {} module on sample: {}\n".format(
                    _format_string_printing_font("clean", ["BOLD"]),
                    _format_string_printing_font(id_key, ["BOLD"]),
                )
            )

            input_bam_path = os.path.join(values[0])
            merged_region = glob.glob(
                os.path.join(os.path.dirname(values[1]))
                + "/{}/*merge_region.tsv".format(id_key)
            )[0]

            cleaned_bam_outdir = os.path.join(values[3])

            _clean(
                script=script,
                bam=input_bam_path,
                sample_id=id_key,
                region=merged_region,
                output_dir=cleaned_bam_outdir,
            )


def _run_make_clean_fastq_batch(bam2fastq_binary, SampleDict, n_cores):

    for [id_key, values] in SampleDict.items():
        cleaned_bam_path = os.path.join(values[3], id_key + ".clean.bam")
        clean_fastq_path = os.path.join(values[2])
        _bam2fastq_10x(
            bam2fastq_binary=bam2fastq_binary,
            bam_file=cleaned_bam_path,
            fastq_outpath=clean_fastq_path,
            n_threads=n_cores,
        )
