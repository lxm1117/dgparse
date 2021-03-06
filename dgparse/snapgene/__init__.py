# encoding=utf-8
"""
Parse Snap Gene File Format and adapt to DTG Schema.
"""
import hashlib
import functools
import uuid

from dgparse import sequtils

from .main import main, parse_snapgene

# Explicitly define our mapping
STRAND = [0, 1, -1, 0]


def extract_sequence(snap_data):
    """
    Adapt a snapgene sequence to a deskgen sequence
    :return:
    """
    bases = snap_data['DNA'].pop('sequence').upper()  # normalize case
    sha1 = hashlib.sha1(bases).hexdigest()
    return {'bases': bases, 'sha1': sha1, 'type_': "dnamoleculesequence"}


def extract_feature_category(snapfeat):
    return snapfeat.pop('type')


def extract_feature(annotation_data, bases):
    """
    Extract the underlying feature
    :return:
    """
    name = annotation_data.pop('name')
    category = extract_feature_category(annotation_data)
    pattern = {'bases': bases, 'sha1': hashlib.sha1(bases).hexdigest()}
    description = annotation_data['Notes'].pop('note', None)
    return {
        'type_': 'dnafeature',
        'name': name,
        'category': category,
        'description': description,
        'sha1': uuid.uuid4().hex,
        'properties': annotation_data['Notes'],
        'pattern': pattern,
        'length': len(pattern['bases']),
    }


def extract_coordinates(annotation_data):
    """
    Split range and return coordinates
    :return:
    """
    # keep segment as it contains additional data
    range_ = annotation_data['Segment']['range'].split('-')
    # may need to subtract one
    directionality = annotation_data.pop('directionality', 0)
    strand = STRAND[directionality]
    start = int(range_[0]) - 1 # convert from [1,n] to pythonic [0,n) coordinate system
    end = int(range_[1])
    return start, end, strand


def extract_annotation(sequence, is_circular, annotation_data):
    """
    Extract the list of annotated features in this sequence
    :return:
    """
    # TODO this object will also need properties to store formatting info
    start, end, strand = extract_coordinates(annotation_data)
    if start < end:
        bases = sequence['bases'][start:end]
    elif is_circular:
        # Feature spans replication origin
        bases = sequence['bases'][start:] + sequence['bases'][:end]
    if not bases:
        raise ParserException('No bases could be parsed for a feature')
    if strand < 0:
        bases = sequtils.get_reverse_complement(bases)
    dnafeature = extract_feature(annotation_data, bases)
    return {
        'type_': 'dnamolecule_dnafeature',
        'start': start,
        'end': end,
        'strand': strand,
        'dnafeature': dnafeature,
    }


def drop_source(annotation_dict):
    dnafeature = annotation_dict.get('dnafeature', None)
    if dnafeature:
        if dnafeature.get('category') != 'source':
            return True


def extract_molecule(molecule):
    """
    Extract the DNA molecule own data
    :param molecule:
    :return:
    """
    is_circular = (molecule['DNA'].pop('topology') == 'circular')
    sequence = extract_sequence(molecule)
    extractor = functools.partial(extract_annotation, sequence, is_circular)
    feature_array = molecule.pop('features', [])
    annotations = map(extractor, feature_array) if feature_array else []
    annotations = filter(drop_source, annotations) 
    description = molecule.pop('descriptor').get('name', None)
    properties = molecule.pop('notes')
    more_props = molecule.pop('other_properties')
    seq_props = molecule.pop('DNA')  # get the remaining items
    properties.update(more_props)
    properties.update(seq_props)
    return {
        '__class__': 'dnamolecule',
        'sequence': sequence,
        'length': len(sequence['bases']),
        'dnafeatures': annotations,  # using updated naming
        'description': description,
        'is_circular': is_circular,
        'properties': properties  # combined properties
    }


def parse(open_file):
    """
    Parse an open snapgene file and convert it into standard DeskGen format
    :param open_file:
    :return:
    """
    snap_result = parse_snapgene(open_file)
    return extract_molecule(snap_result)

