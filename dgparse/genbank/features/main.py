import logging
import re
import string

from .constants import FEATURE_TYPES, FEATURE_QUALIFIERS
from ..constants import GENBANK_HEADERS
from dgparse.exc import NullCoordinates
logger = logging.getLogger(__name__)


def parse_coord(coord):
    """Parse the coordinates"""
    # For full syntax see http://www.ddbj.nig.ac.jp/FT/full_index.html#3.4
    # Suffice here to parse the max and min values
    coord_dict = {'strand': -1 if "complement" in coord else 1}    
    start = None
    end = None
    for pos in re.findall(r'\d+', coord):
        if start is None or int(pos) < start:
            start = int(pos)
        if end is None or int(pos) > end:
            end = int(pos)
    start -= 1 # convert from [1,n] to pythonic [0,n) coordinate system
    coord_dict.update({'start': start, 'end': end})
    if (end - start) < 1:
        raise NullCoordinates
    return coord_dict


def recurse_features(line, lines, out):
    'Handle blank lines, manage feature parsing'
    tokens = line.strip().split()
    if tokens[0] in GENBANK_HEADERS:
        return line, lines, out
    # Offsets per http://www.ddbj.nig.ac.jp/FT/full_index.html#3.4
    feature_key = line[5:20].strip()
    feature_content = line[21:].strip()
    qualifier_match = re.match(r'/(\w+)=\"(.*)', feature_content)
    if feature_key:
        # New feature, instantiate, initialize qualifier, append
        out['features'].append({'category': feature_key})
        out['features'][-1].update(parse_coord(feature_content))
        if feature_key not in FEATURE_TYPES:
            logger.info("Feature '{0}' not recognised".format(feature_key))
    elif qualifier_match:
        # New qualifier, initialize and append
        qualifier_key = qualifier_match.group(1)
        qualifier_val = qualifier_match.group(2).rstrip('"')
        out['features'][-1].update({qualifier_key: qualifier_val})
        if qualifier_key not in FEATURE_QUALIFIERS:
            logger.info("Qualifier '{0}' not recognised".format(qualifier_key))
        # Assign first qualifier as default feature name
        if 'name' not in out['features'][-1]:
            out['features'][-1].update({'name': qualifier_val})
    else:
        pass
        # Multiline qualifier, ignore remainder
    return recurse_features(next(lines), lines, out)
