
import hashlib

from .main import main, parse_snapgene
from ..exc import ParserException



def parse(open_file):
    """
    Parse an open snapgene file and convert it into standard DeskGen format
    :param open_file:
    :return:
    """
    result = parse_snapgene(open_file)
    try:
        snap_dna = result.pop('DNA')
        bases = snap_dna.pop('sequence').upper()  # normalize case
        properties = result.pop('other_properties')
        properties.update()
        # adapt topology
        result['is_circular'] = (snap_dna.get('topology', 'circular')
                                 == 'circular')

        # remove 'source' as these are not real features.
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
        },
        'properties': properties,
    })
    return result
