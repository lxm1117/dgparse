import hashlib
from . import main
from ..exc import ParserException


def parse(open_file):
    'Schema compatibility for new parser'
    result = main.init(open_file)
    try:
        bases = result.pop('origin')
        def drop_source(feature_dict):
            if feature_dict.get('category') != 'source':
                return True
        result['features'] = filter(drop_source, result['features'])
    except KeyError:
        raise ParserException('No sequence could be parsed')
    result.update({
        'sequence': {
            'bases': bases,
            'sha1': hashlib.sha1(bases).hexdigest(),
        }
    })
    return result
