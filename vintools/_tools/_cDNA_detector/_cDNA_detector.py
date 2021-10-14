# package imports #
# --------------- #
import os

# local imports #
# ------------- #
from ._supporting_functions._tabulate_results import _tabulate_results
from ._supporting_functions._make_10x_SampleDict import (
    _organize_10x_samples,
    _add_cleaned_outpath_to_10x_SampleDict,
    _find_orphaned_samples,
)
from ._supporting_functions._view import _view
from ._supporting_functions._main_funcs import (
    _detect,
    _format_outdir,
    _prepend_gene_model,
    _prepare,
    _setup_cDNA_detector,
    _bam2fastq_10x,
    _run_make_clean_fastq_batch,
    _download_bam2fastq_10x_binary,
    _clean_bam_batch,
    _run_clean_preflight,
)

from ..._utilities._system_utils._fetch_cpu_count import _fetch_cpu_count
from ..._utilities._ux_utils._pystrings import _format_string_printing_font
from ..._utilities._system_utils._get_pypi_package_loc import _get_pypi_package_loc


class _cDNA_detector:
    def __init__(self, script=False, n_cores=False, bam2fastq_binary=False):

        """
        Instantiate cDNA Detector. 
        
        Parameters:
        -----------
        script
            type: str or bool
            
        n_cores
            type: int or bool
        
        Returns:
        --------
        None
        
        Notes:
        ------
        (1)  By default, all cores are used, unless otherwise specified with `n_cores`.
        (2)  If the GitHub repository is not found, it is cloned and installed with its respective
             dependencies.
        (3)  The `__init__()` function of this class clones the cDNA-detector repository adjacent to the 
             repository from which I am executing the commands in this wrapper. Thus, this wrapper references
             the main python script in that library: `/home/mvinyard/software/cDNA-detector/cdna-detector.py`
        """

        if not script:
            software_dir = os.path.dirname(_get_pypi_package_loc())
            cDNA_detector_repo_path = os.path.join(software_dir, "cDNA-detector")
            self.script = os.path.join(cDNA_detector_repo_path, "cdna-detector.py")
            if not os.path.exists(self.script):
                _setup_cDNA_detector(cDNA_detector_repo_path, script)
        if not bam2fastq_binary:
            self.bam2fastq_binary = os.path.join(
                os.path.dirname(_get_pypi_package_loc()), "bamtofastq_linux"
            )

        if not n_cores:
            self.n_cores = _fetch_cpu_count()

        self.clean_preflight_complete = False

    def prepare(self, ref_10x):

        """
        Prepare a gene model from a 10x reference. 
        
        Parameters:
        -----------
        ref_10x
        
        Returns:
        --------
        Prepends gene model to class. 
        
        Notes:
        ------
        """

        self.gene_model = _prepare(self.script, ref_10x)

    def get_samples(self, data_path=False, outdir=False):

        """
        """

        if data_path:
            self.data_path = data_path
        if outdir:
            self.outdir = outdir

        self.SampleDict = _organize_10x_samples(self.data_path, self.outdir)

    def preflight(self, data_path=False, gene_model=False, outdir=False):

        """
        Run the preflight module. 
        
        Parameters:
        -----------
        data_path
            type: str
        
        gene_model
            type: str or bool
        
        outdir
            type: str or bool
        
        Returns:
        --------
        None, defines the following:
            
            self.data_path
            self.gene_model
            self.SampleDict
        
        Notes:
        ------
        
        """

        if not data_path:
            pass
        else:
            self.data_path = data_path

        _format_outdir(self, outdir)
        if gene_model:
            self.gene_model = gene_model
        else:
            print(
                "Using available gene model:\n\n\t{}".format(
                    _format_string_printing_font(self.gene_model, ["BOLD", "BLUE"])
                )
            )

        self.SampleDict = _organize_10x_samples(self.data_path, self.outdir)

    def view(self, orphaned=False, data_path=False, outdir=False):

        """
        prints self.SampleDict contents in readable format. 
        
        Parameters:
        -----------
        None, self
        
        data_path [ optional ]
        
        outdir [ optional ]
        
        Returns:
        --------
        prints self.SampleDict contents in readable format. 
        
        Notes:
        ------
        
        """

        try:
            self.SampleDict
        except:
            if data_path:
                self.data_path = data_path
            if outdir:
                self.outdir = outdir
            self.SampleDict = _organize_10x_samples(self.data_path, self.outdir)

        if not orphaned:
            _view(self.SampleDict)
        else:
            _view(self.OrphanedSampleDict)

    def detect(
        self, orphaned=False, fast=False, dry_run=False,
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

        if not orphaned:
            Dict = self.SampleDict
        else:
            Dict = self.OrphanedSampleDict

        for [sample, values] in Dict.items():
            _detect(
                script=self.script,
                bam=values[0],
                sample_id=sample,
                gene_model=self.gene_model,
                n_cores=self.n_cores,
                fast=fast,
                cDNA_detector_outdir_sample=values[1],
                dry_run=dry_run,
            )

    def get_orphaned_samples(self):

        """
        
        """

        self.OrphanedSampleDict = _find_orphaned_samples(self.SampleDict)

    def tabulate_results(self):

        """
        Tabulate and plot results.
        
        Parameters:
        -----------
        None, self
        
        Returns:
        --------
        None, defines: 
        
            self.results
        """

        self.ResultsDict, self.missing_results = _tabulate_results(self.outdir)

    def clean_preflight(self):

        """
        
       
        """

        _run_clean_preflight(self)

    def clean(self):

        """
        
        """

        if self.clean_preflight_complete == False:
            _run_clean_preflight(self)
        _clean_bam_batch(self.script, self.SampleDict, self.OrphanedSampleDict)

    def make_fastq(self):

        """"""

        if not os.path.exists(self.bam2fastq_binary):
            self.bam2fastq_binary = _download_bam2fastq_10x_binary()

        _run_make_clean_fastq_batch(
            bam2fastq_binary=self.bam2fastq_binary,
            SampleDict=self.SampleDict,
            n_cores=self.n_cores,
        )
