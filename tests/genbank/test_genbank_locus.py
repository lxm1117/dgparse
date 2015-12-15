from dgparse.genbank.locus import parse as parse_locus

def test_parsers_genbank_locus(genbank_locus_lines):
    'Genbank locus lines are happily parsed into dicts'
    input_output = zip(*genbank_locus_lines)
    for input, expected_output in input_output:
        output = parse_locus(input, {})['locus']
        tuple_output = (
            output['name'],
            output['length'],
            output['molecule_type'],
            output['division'],
            output['orientation'],
            output['accessed'],
        )
        assert expected_output == tuple_output
