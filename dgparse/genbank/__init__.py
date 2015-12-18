import hashlib
import copy
from . import main
from ..exc import ParserException

from dgparse import sequtils

def pick_a_name(dict_):
    """Pick a name from a dnafeature dict. The rules are as follows:
    1) Label or Name ==> name
    2) The first occurence of note ==> label ==> name
    3) "Novel" as in a new feature of a known type.
    """
    # label, then note, then name, then default
    for key in ['label', 'note', 'name']:
        if key in dict_:
            return dict_.pop(key)
    return 'Novel'


def pick_description(dict_):
    """Pick a description from a dna feature dict from the parser.
    1) Note after label
    2) The second occurence of "note" for a feature
    """
    default = "Automatically imported"
    try:
        note = dict_.pop('note') + '\n'
    except KeyError:
        note = ''
    return note + default


def drop_source(feature_dict):
    if feature_dict.get('category') != 'source':
        return True


def parse(open_file):
    'Parse an open genbank file and convert it into standard DeskGen format'
    result = main.init(open_file)
    try:
        bases = result.pop('origin')
    except KeyError:
        raise ParserException('No sequence could be parsed')
    result.update({
        'sequence': {
            'bases': bases,
            'sha1': hashlib.sha1(bases).hexdigest(),
        }
    })
    try:
        result['is_circular'] = result['locus']['orientation']
    except KeyError:
        raise ParserException('No orientation could be parsed')
    result['dnafeatures'] = list() # a list of dnafeature annotations
    features = result.get('features', {}) # TODO: pop this, replaced by 'dnafeatures'
    features = filter(drop_source, features) 
    for feature in features:
        print "Feature %s"%(feature)
        unpack = copy.deepcopy(feature)
        annotation = dict()
        for key in 'start', 'end', 'strand':
            annotation[key] = unpack.pop(key)
        if annotation['start'] < annotation['end']:
            bases = result['sequence']['bases'][annotation['start']:annotation['end']]
        elif result['is_circular']:
            # Feature spans replication origin
            bases = result['sequence']['bases'][annotation['start']:] + \
                    result['sequence']['bases'][:annotation['end']]
        if annotation['strand'] < 0:
            bases = sequtils.get_reverse_complement(bases) # assumed pythonic coordinates
        annotation['dnafeature'] = {
            'name': pick_a_name(unpack),
            'category': unpack.pop('category', None),
            'description': pick_description(unpack),
            'length': len(bases),
            'pattern': {
                'bases': bases,
                'sha1': hashlib.sha1(bases).hexdigest()
            },
        }
        annotation['dnafeature']['properties'] = unpack # anything else
        result['dnafeatures'].append(annotation)
    return result
