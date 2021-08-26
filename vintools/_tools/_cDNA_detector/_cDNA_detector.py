
# package imports #
# --------------- #
import os

# local imports #
# ------------- #
from ._supporting_functions._tabulate_results import _tabulate_results
from ._supporting_functions._make_10x_SampleDict import _organize_10x_samples
from ._supporting_functions._main_funcs import _detect, _view, _format_outdir, _prepend_gene_model, _prepare

from ..._utilities._fetch_cpu_count import _fetch_cpu_count
from ..._utilities._fetch_cpu_count import _fetch_cpu_count
from ..._utilities._pystrings import _format_string_printing_font
from ..._utilities._clone_GitHub_repo import _clone_GitHub_repo

class _cDNA_detector:
    def __init__(self, script=False, n_cores=False):

        """"""
        
        print(os.getcwd())
        print(self.__path__)

#         if not script:
#             self.script = "../../../../cDNA-detector/cdna-detector.py" # navigates out of vintools into the adjacent repo
#             if not os.path.exists(self.script):
#                 print("cDNA-detector repository not found... cloning locally now.")
#                 _clone_GitHub_repo(repository="https://github.com/rheinbaylab/cDNA-detector.git",
#                                    destination="../../../../cDNA-detector")
#         if not n_cores:
#             self.n_cores = _fetch_cpu_count()
            
            
    def prepare(self):
    
        """"""
    
        self.gene_model = _prepend_gene_model(OUTDIR)
        _prepare(script, ref_10x, outdir=os.path.dirname(self.gene_model))

    def detect_preflight(self, data_path, gene_model, outdir=False):

        """"""

        self.data_path = data_path
        self.gene_model = gene_model
        _format_outdir(self, outdir)

        self.SampleDict = _organize_10x_samples(data_path, self.outdir)
        self.gene_model = gene_model
    
    def view(self):
        
        """
        prints self.SampleDict contents in readable format. 
        
        Parameters:
        -----------
        None, self
        
        Returns:
        --------
        prints self.SampleDict contents in readable format. 
        
        Notes:
        ------
        
        """
        
        _view(self)
        
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
            
    def tabulate_results(self):

        """Tabulate and plot results."""

        self.results = _tabulate_results(self.outdir)