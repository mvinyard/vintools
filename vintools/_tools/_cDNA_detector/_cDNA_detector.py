
# local imports #
# ------------- #
from ._supporting_functions._make_10x_SampleDict import _organize_10x_samples
from ._supporting_functions._main_funcs import _detect, _view, _format_outdir

from ..._utilities._fetch_cpu_count import _fetch_cpu_count
from ..._utilities._fetch_cpu_count import _fetch_cpu_count
from ..._utilities._pystrings import _format_string_printing_font


class _cDNA_detector:
    def __init__(self, script=False, n_cores=False):

        """"""

        if not script:
            self.script = "/home/mvinyard/software/cDNA-detector/cdna-detector.py"
        if not n_cores:
            self.n_cores = _fetch_cpu_count()

    def preflight(self, data_path, gene_model, outdir=False):

        """"""

        self.data_path = data_path
        self.gene_model = gene_model
        _format_outdir(self, outdir)

        self.SampleDict = _organize_10x_samples(data_path, self.outdir)
        self.gene_model = gene_model
    
    def view(self):
        
        _view(self)
        
        for [sample, [bam_in, out]] in self.SampleDict.items():
            
            print("Sample: {}\n{bam_in}\n{out}\n".format(_format_string_printing_font(sample, ['BOLD'])))
            
        
    def detect(
        self,
        dry_run=False,
    ):

        """
        Main detect module of cDNA-detector. 
        
        Parameters:
        -----------
        dry_run
            default: False
            type: bool
            
        Returns:
        --------
        None
            outputs to designated directory
            
        Notes:
        ------
        (1) Must instantiate cDNA-detector class and run preflight prior to running this function.
        
        For more, see adjacent README.md
        """

        for [sample, [bam_in, out]] in self.SampleDict.items():
            _detect(
                script=self.script,
                bam=bam_in,
                sample_id=sample,
                gene_model=self.gene_model,
                n_cores=self.n_cores,
                cDNA_detector_outdir_sample=out,
                dry_run=dry_run,
            )