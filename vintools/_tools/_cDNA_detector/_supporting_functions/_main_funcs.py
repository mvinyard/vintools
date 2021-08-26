
# package imports #
# --------------- #
import os
import glob

# local imports #
# ------------- #
from ...._utilities._pystrings import _format_string_printing_font


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

def _view(self):

    for [sample, [bam_in, out]] in self.SampleDict.items():

        print(
            "{}: {}\n\n\t{}\n\t{}\n".format(
                _format_string_printing_font("Sample", ["BOLD"]),
                _format_string_printing_font(sample, ["BOLD", "RED"]),
                _format_string_printing_font(bam_in, ["BOLD"]),
                _format_string_printing_font(out, ["BOLD"]),
            )
        )

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
        
def _prepend_gene_model(prepare_outdir):
    
    """"""
    
    gene_model = glob.glob(prepare_outdir + "*.saf")[0]
    
    return gene_model

def _prepare(script, ref_10x, outdir):
    
    """"""        

    annotations=os.path.join(ref_10x, "genes/genes.gtf")
    genome=os.path.join(ref_10x_10x, "fasta/genome.fa")
    
    executable = "{} prepare --annotation {} --genome {} --output_dir {}".format(script, annotations, genome, outdir)
    os.system(executable)
    
    _prepend_gene_model(prepare_outdir)