import os
from dgparse.fasta.parse_fasta import parse_fasta_str


def test_basic_fasta_file():
    path = os.path.join(os.path.dirname(__file__), "../../data/fasta/ppiczalphaa.fasta")
    with open(path, 'r') as test_file:
        test_str = test_file.read()
        output = parse_fasta_str(test_str)
        assert output
