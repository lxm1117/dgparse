import logging
import re
import string

from .constants import FEATURE_TYPES, FEATURE_QUALIFIERS
from ..constants import GENBANK_HEADERS
from dgparse.exc import NullCoordinates
logger = logging.getLogger(__name__)


def parse_coord(coord):
    """Parse the coordinates"""
    coord_dict = {'strand': -1 if "complement" in coord else 1}
    start, end = map(int, coord.strip('complement()\r\n').split('..')) # !! operator syntax is more extensive see http://www.ddbj.nig.ac.jp/FT/full_index.html#3.4
    start -= 1 # convert from [1,n] to pythonic [0,n) coordinate system
    coord_dict.update({'start': start, 'end': end})
    if (end - start) < 1:
        raise NullCoordinates
    return coord_dict


def recurse_qualifier(line, lines, feature, qualifier):
    break_conditions = (
        any(map(lambda H: H in line, GENBANK_HEADERS)),  # next header
        any(map(lambda F: F in line, FEATURE_TYPES)),  # next feature
    )
    if any(break_conditions):
        return line, lines, feature
    if bool(re.match(r'^/\w+=', line.strip())):  # next qualifier
        return recurse_feature(line, lines, feature)
    else:
        feature[qualifier] = ' '.join([feature[qualifier],
                                       line.strip(' "\r\n')])
    return recurse_qualifier(next(lines), lines, feature, qualifier)


def recurse_feature(line, lines, feature):
    'Handle blank lines in features, manage qualifier parsing'
    tokens = line.strip().split()
    if not tokens:
        return recurse_feature(next(lines), lines, feature)
    qualifier_match = re.match(r'^(/\w+=)', line.strip())
    if qualifier_match:
        qualifier_key = qualifier_match.group(0)
        if qualifier_key not in FEATURE_QUALIFIERS:
            logger.info("Qualifier '{0}' not recognised".format(qualifier_key))
        qualifier = qualifier_key.strip('/=')
        qualifier_value = line.split('=')[-1].strip('"\r\n')
        # Don't lose info on the double-note case
        if qualifier == 'note' and 'note' in feature:
            # rename the first note to the label, so it get's mapped to name
            feature['label'] = feature['note']
            # the new note takes the 'note' attribute
            feature['note'] = qualifier_value
        else:
            feature.update({qualifier: qualifier_value})
        return recurse_qualifier(next(lines), lines, feature, qualifier)
    return line, lines, feature


def recurse_features(line, lines, out):
    'Handle blank lines, manage feature parsing'
    tokens = line.strip().split()
    if not tokens:
        return recurse_features(next(lines), lines, out)
    if tokens[0] in GENBANK_HEADERS:
        return line, lines, out
    if len(tokens) == 1:  # Throw it away
        return recurse_features(next(lines), lines, out)
    if tokens[0].lower() in map(string.lower, FEATURE_TYPES):
        feature = {'category': tokens[0]}
        try:
            feature.update(parse_coord(tokens[1]))
        except:
            return recurse_features(next(lines), lines, out)
        line, lines, feature = recurse_feature(next(lines), lines, feature)
        out['features'].append(feature)
        tokens = line.strip().split()
        if tokens[0] in GENBANK_HEADERS:  # fff
            return recurse_features(line, lines, out)
    return recurse_features(line, lines, out)
