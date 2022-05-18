import pytest
from Bio import AlignIO
from probedesign import probesearch

def test_get_sequence_regions():
    test_alignment = probesearch.alignment("tests/test_align.fasta")
    sequence_regions = test_alignment._get_sequence_regions()
    assert(sequence_regions["Pratylenchus-coffeae"] == (0, 1611))