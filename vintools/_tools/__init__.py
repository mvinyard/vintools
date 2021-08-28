
# tools __init__.py

__author__ = ', '.join([
    'Michael E. Vinyard'
])
__email__ = ', '.join([
    'vinyard@g.harvard.edu',
])


from ._cDNA_detector._cDNA_detector import _cDNA_detector as cDNA_detector

from ._csv_to_markdown_table._csv_to_markdown_table import create_md_table
from ._csv_to_markdown_table._csv_to_markdown_table import quicktable

from ._DESeq2._batch_DESeq2 import _run_batch_DESeq2 as DESeq2
