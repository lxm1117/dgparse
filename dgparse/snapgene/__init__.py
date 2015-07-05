
import hashlib
from .main import main, parse_snapgene
from ..exc import ParserException

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
    feature_type = snapfeat.pop('type')
    return {
        'type_': 'dnafeaturecategory',
        'name': feature_type,
        'accession': feature_type.lower(),  # normalize caps
    }


def extract_feature(annotation_data):
    """
    Extract the underlying feature
    :return:
    """
    name = annotation_data.pop('name')
    category = extract_feature_category(annotation_data)
    if category['accession'] == 'source':
        return  # No source features
    accession = category['name'] + '-' + name
    return {
        'type_': 'dnafeature',
        'name': name,
        'category': category,
        'accession': accession,
        'properties': annotation_data
    }

def extract_coordinates(annotation_data):
    """
    Split range and return coordinates
    :return:
    """
    # keep segment as it contains additional data
    range = annotation_data['Segment']['range'].split('-')
    # may need to subtract one
    directionality = annotation_data.pop('directionality', 0)
    start = int(range[0])
    end = int(range[0])
    return start, end, STRAND[directionality]

def extract_annotation(annotation_data):
    """
    Extract the list of annotated features in this sequence
    :return:
    """
    # TODO this object will also need properties to store formatting info
    start, end, strand = extract_coordinates(annotation_data)
    return {
        'type_': 'dnamolecule_dnafeature',
        'start': start,
        'end': end,
        'strand': strand,
        'dnafeature': extract_feature(annotation_data),
    }

def extract_molecule(molecule):
    """
    Extract the DNA molecule own data
    :param molecule:
    :return:
    """
    is_circular = (molecule['DNA'].pop('topology') == 'circular')
    annotations = map(extract_annotation, molecule.pop('features'))
    description = molecule.pop('descriptor')
    properties = molecule.pop('notes')
    more_props = molecule.pop('other_properties')
    sequence = extract_sequence(molecule)
    seq_props = molecule.pop('DNA')  # get the remaining items
    properties.update(more_props)
    properties.update(seq_props)
    return {
        'type_': 'dnamolecule',
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

