
from ._generate_10x_cellranger_ResultsDict import _generate_10x_cellranger_ResultsDict
from ._compile_web_summary_files import _cp_web_summaries

class _post_10x_pipestance:
    
    """
    Module for easy wrap-up after a 10x cellranger run. Should work with most versions/modalities of cellranger.
    """
    
    def __init__(self):
        
        """"""
        
        
    def make_ResultsDict(self, path_to_results, samples):
        
        """
        Make a dictionary of the results. 
        
        Parameters:
        -----------
        path_to_results
        
        samples
        
        Returns:
        --------
        None
        
        Notes:
        ------
        
        """
        
        self.path_to_results = path_to_results
        self.samples = samples
        self.ResultsDict = _generate_10x_cellranger_ResultsDict(path_to_results, samples)
        
    def cp_web_summaries(self, gcp_bucket, title=False):
        
        """
        Use gsutils to copy web summaries from all samples into a combined folder for storage and sharing. 
        
        Definition:
        -----------
        gcp_bucket
            Google cloud bucket. 
        
        title
            Title of the compiled web summary directory moved a google bucket
            
        
        Returns:
        --------
        None
        
        Notes:
        ------
        (1) One must already be logged in using `gcloud auth login`
        """
        
        _cp_web_summaries(self.ResultsDict, gcp_bucket, title)
        
        
def _run_10x_post_pipestance(path_to_results, samples, gcp_bucket, compiled_web_summary_title=False):
    
    """
    Function for easy wrap-up after a 10x cellranger run. Should work with most versions/modalities of cellranger.
    
    Parameters:
    -----------
    
    Returns:
    --------
    
    Notes:
    ------
    
    """
    
    post10x = post_10x_pipestance()
    
    post10x.make_ResultsDict(path_to_results, samples)
    post10x.cp_web_summaries(gcp_bucket, compiled_web_summary_title)