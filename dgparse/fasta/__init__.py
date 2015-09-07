from . import parse_fasta
from dgparse.exc import ParserException


def parse(open_file):
    'Interface compatibility for old fasta parser'
    result = parse_fasta.parse_fasta_file(open_file)
    if result:
        result = result[0]
        result.update({
            'sequence': {
                'bases': result['sequence']['seq'],
                'sha1': result['sequence']['sha1'],
            }
        })
        return result
    raise ParserException('Fasta parse fail!')
