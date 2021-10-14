# _genome_editing __init__.py

__author__ = ", ".join(["Michael E. Vinyard"])
__email__ = ", ".join(["vinyard@g.harvard.edu",])

# sequence manipulation functions
from ._sequence_manipulation._get_reverse_complement import (
    _get_reverse_complement as reverse_complement,
)
from ._sequence_manipulation._get_chromosome_sequence import (
    _get_chromosome_sequence as get_chromosome_sequence,
)
