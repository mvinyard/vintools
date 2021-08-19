
# local imports #
# ------------- #
from ._cDNA_detector_detect import _detect, _organize_10x_samples, _format_outdir
from ..._utilities._fetch_cpu_count import _fetch_cpu_count

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

    def detect(
        self,
        dry_run=False,
    ):

        """
        Parameters:
        -----------
        
        Returns:
        --------
        
        Notes:
        ------
        
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