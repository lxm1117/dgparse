import tempfile
from . import locus
from .features import recurse_features
from .origin import recurse_origin


def recurse_headers(line, lines, out):
    '''
    Try to select an appropriate function to continue parsing with, otherwise
    skip the line.
    '''
    tokens = line.strip().split()
    if not tokens:
        return recurse_headers(next(lines), lines, out)
    func = HEADER_FUNCTIONS.get(tokens[0], False)
    if func is False:
        return recurse_headers(next(lines), lines, out)
    return func(line, lines, out)


def parse_features(line, lines, out):
    'Parse features from lines, appending to out dict'
    out['features'] = list()
    line, lines, out = recurse_features(next(lines), lines, out)
    return recurse_headers(line, lines, out)


def parse_origin(line, lines, out):
    out['origin'] = ''
    line, lines, out = recurse_origin(next(lines), lines, out)
    return recurse_headers(line, lines, out)


def parse_locus(line, lines, out):
    'Call locus parser, return iterator to main thread'
    out = locus.parse(line, out)
    return recurse_headers(next(lines), lines, out)


HEADER_FUNCTIONS = {
    'LOCUS': parse_locus,
    #'DEFINITION'
    #'ACCESSION'
    #'VERSION'
    #'KEYWORDS'
    #'SOURCE'
        #'ORGANISM'
    #'REFERENCE'
        #'AUTHORS'
        #'TITLE'
        #'JOURNAL'
        #'PUBMED'
    'FEATURES': parse_features,
    'ORIGIN': parse_origin,
}


def safe_file(open_file):
    'Reopen file with U for universal newlines, return a generator'
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(open_file.read())
        with open(tmp.name, 'U') as safe_tmp:
            for line in safe_tmp:
                yield line


def init(open_file, default=None):
    'Parse an open GenBank file to a dict, init for recursion'
    # lines = safe_file(open_file)
    lines = iter(open_file.readlines())
    out = dict()
    try:
        # Trust good old call-by-reference to do the job
        recurse_headers(next(lines), lines, out)
    except StopIteration:
        return out
