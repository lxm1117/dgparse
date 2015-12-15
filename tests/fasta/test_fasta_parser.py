'''Test the genbank parser works correctly for several cases'''
import os
import pytest
from dgparse.fasta.parse_fasta import parse_fasta_str


@pytest.mark.parametrize('path',(
    '../../data/fasta/ppiczalphaa.fasta',
    '../../data/fasta/Cellecta-pRSI16-U6-sh-HTS6-UbiC-TagRFP-2A-Puro.fasta',
))
def test_parse_fasta_str(path):
    """
    GIVEN an input file and a file format sequence, return the correct
    sequence record and sequence object.
    """
    abpath = os.path.join(os.path.dirname(__file__), path)
    result = parse_fasta_str(open(abpath, 'r').read())
    assert len(result) > 0
