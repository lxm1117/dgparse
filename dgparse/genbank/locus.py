# encoding=utf-8
"""
Parse the locus line of a genbank file
"""
from . import constants
from datetime import datetime
from dgparse.exc import ParserException


def parse_known_date_formats(potential_date_str):
    'Parse a string by known date formats, return the first one or None.'
    def try_date_format(date_format):
        try:
            return datetime.strptime(potential_date_str, date_format)
        except ValueError:
            return False
    results = filter(None, map(try_date_format,
                               constants.KNOWN_DATE_FORMATS))
    if any(results):
        return results[0]


def look_back(tokens, div=None, ori=None, expect=False):
    if tokens[-1] in constants.GENBANK_DIVISONS:
        div = tokens.pop(-1)
    elif tokens[-1] in constants.VALID_ORIENTATIONS:
        ori = tokens.pop(-1)
    if not any([div, ori]) and expect:
        message = ('No orientation or division found in LOCUS line:'
                   "\n{0}".format(' '.join(tokens)))
        raise ParserException(message)
    return tokens, div, ori


def division_or_orientation_or_both(header_tokens, three_chances):
    """
    When passed list of header token, get either a Genbank division or an
    orientation or both.
    Return the list with consumed items removed.
    """
    # look at the first token, we don't expect a match here if there was a
    # malformed date in the previous frame
    header_tokens, division, orientation = look_back(header_tokens,
                                                     expect=not three_chances)
    if not any([division, orientation]) and three_chances:
        header_tokens.pop(-1)
    # look at the second token
    header_tokens, division, orientation = look_back(header_tokens,
                                                     division,
                                                     orientation)
    # if we failed to find a date, the third token back could be div/ori
    if three_chances:
        # look at the third token
        header_tokens, division, orientation = look_back(header_tokens,
                                                         division,
                                                         orientation)
    return (
        header_tokens,
        division,
        constants.ORIENTATION_MAPPING.get(orientation),
    )


def parse(line, out):
    """Have a good go at parsing data from the LOCUS line"""
    locus_tokens = line.split()
    potential_date_str = locus_tokens[-1]
    accessed = parse_known_date_formats(potential_date_str)
    if accessed:
        locus_tokens.pop(-1)
    (locus_tokens,
     division,
     orientation) = division_or_orientation_or_both(locus_tokens,
                                                    not bool(accessed))
    molecule_type = locus_tokens.pop(-1)
    bp = locus_tokens.pop(-1)
    if bp == 'bp':
        length = int(locus_tokens.pop(-1))
    else:
        length = None
    if locus_tokens.pop(0) != 'LOCUS':
        msg = 'LOCUS was not at the start of the LOCUS line in file'
        raise ParserException(msg)
    name = ' '.join(locus_tokens)
    out['locus'] = {
        'name': name,
        'length': length,
        'molecule_type': molecule_type,
        'division': division,
        'orientation': orientation,
        'accessed': accessed,
    }
    return out
