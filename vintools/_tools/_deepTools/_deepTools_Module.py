

# package imports #
# --------------- #
import os

# local imports #
# ------------- #
from ..._utilities._system_utils._use_n_cores import _use_n_cores
from ..._utilities._system_utils._flexible_mkdir import _flexible_mkdir
from ._deepTools_supporting_functions import _find_sample_files, _plotCorrelation

class _deepTools:
    
    def __init__(self, silent=False, verbose=False):
        
        """Notebook interface of deepTools"""
        
        self.silent = silent
        self.verbose = verbose
        self.DataDict = {}
        self.DataDict['bams'] = {}
        self.DataDict['bigwigs'] = {}
        
    def load_bams(self, path, sample_dirname=False, silent=False):

        if silent:
            self.silent = silent
        
        self.DataDict['bams'].update(_find_sample_files(path, file_extension="bam", sample_dir_level=sample_dirname, silent=silent))
    
    def load_bw(self, path, sample_dirname=False, silent=False):

        if silent:
            self.silent = silent
        
        self.DataDict['bigwigs'].update( _find_sample_files(path, file_extension="bw", sample_dir_level=sample_dirname, silent=silent))
    
    def bamCoverage(self, BamDict=False, output_dir="./", dry_run=False, verbose=False, n_cores=2):
        
        """If supplied, a dictionary of format: {'sample': '/path/to/sample/sample.bam'} is taken as input"""
        
        if verbose:
            self.verbose = verbose
        
        self.n_cores = _use_n_cores(n_cores)
        
        print("\tnote... this is {} cores per sample. Samples are run as background processes for parallelization.\n".format(self.n_cores))
        
        if BamDict:
            self.BamDict = BamDict
            
        if not os.path.exists(output_dir):
            _flexible_mkdir(output_dir)
            
        if dry_run:
            print("deepTools bamCoverage command issued:\n")
            
        for sample, bamfile in self.BamDict.items():
            outprefix = os.path.join(output_dir, sample)
            self.bamCoverage_executable = "bamCoverage -b {} -o {}.bw -p {} -of bigwig &".format(bamfile, outprefix, n_cores)
            if dry_run:
                print(self.bamCoverage_executable, "\n")
            else:
                if self.verbose:
                    print("deepTools bamCoverage command issued:", "\n")
                    print(self.bamCoverage_executable, "\n")
                os.system(self.bamCoverage_executable)
                
    def multiBigwigsummary(self, outfile="multiBigwigsummary.results.npz", n_cores=False):
    
        """runs multiBigwigsummary from deepTools."""
        
        self.outfile = outfile
        self.outdir = os.path.dirname(self.outfile)
        if not os.path.exists(self.outdir):
            _flexible_mkdir(self.outdir)
            

        bigwigs = " ".join(list(self.DataDict['bigwigs'].values()))
        
        self.n_cores = _use_n_cores(n_cores)
        
        self.multiBigwigsummary_executable = "multiBigwigSummary bins -b {} -o {} -p {}".format(bigwigs, outfile, self.n_cores)
        os.system(self.multiBigwigsummary_executable)
    
    def plotCorrelation(self, title=False, summary_results_file=False, silent=False):
        
        if summary_results_file:
            self.outfile = outfile
            self.outdir = os.path.dirname(self.outfile)
               
        if not title:
            title=self.outdir
        else:
            title=os.path.join(self.outdir, title)
                
        _plotCorrelation(title, self.outfile, silent)