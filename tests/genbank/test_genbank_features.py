# -*- coding: utf-8 -*-
"""
Test genbank feature annotations
"""

import os
import io
import csv
import pytest

from dgparse import genbank

@pytest.mark.parametrize("file_name,expected_feature_count", [
    ('01-lentiCRISPRv2-add.gb', 19),
    ('02-pIB2-SEC13-mEGFP-snap.gb', 8),
    ('03-NM_024674-add-construct.gb', 22),
    ('04-px330-snap.gb', 13),
    ('05-lentiCas9-Blast-add.gb', 15),
])
def test_parse_genbank_multifeature_count(file_name, expected_feature_count):
    """Obtain expected feature count for some example files"""
    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/' + file_name)
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        parsed_annotations = ret['dnafeatures']
        print "Names parsed"
        for anno in parsed_annotations:
            print anno['dnafeature']['category'], anno['dnafeature']['name']
        assert len(parsed_annotations) == expected_feature_count
    

def test_parse_genbank_multifeature():
    """Can parse genbank file and retrieve expected features"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/PX330.gbk')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        parsed_annotations = ret['dnafeatures']

    expected_features_tsv = '''\
category	start	end	strand	length	name
misc_feature	439	1238	1	799	CBh
misc_feature	1370	5471	1	4101	hSpCsn1
misc_feature	6311	6543	1	232	bGH polyA
misc_feature	6551	6692	1	141	R-ITR
rep_origin	6783	7090	1	307	f1 Origin
misc_feature	7608	8466	1	858	AmpR
rep_origin	8616	9284	1	668	pUC ori
misc_feature	8616	9284	-1	668	made up rev pUC ori
promoter	0	249	1	249	Human U6 Promoter
misc_feature	5471	5519	1	48	NLS
misc_feature	267	343	1	76	chimeric guide RNA scaffold
misc_feature	343	349	1	6	U6 terminator
misc_feature	1250	1319	1	69	3xFLAG
misc_feature	1319	1370	1	51	NLS
misc_feature	5519	5525	1	6	EcoRI
misc_feature	5525	5588	1	63	T2A
misc_feature	5588	6302	1	714	GFP
misc_feature	6302	6308	1	6	EcoRI
'''
    stream = io.BytesIO()
    stream.write(expected_features_tsv)
    stream.seek(0)
    reader = csv.DictReader(stream, delimiter='\t')
    expected_annotations = list()
    for annotation_dict in reader:
        expected_annotations.append(annotation_dict)

    assert len(parsed_annotations) == len(expected_annotations)

    for expected, parsed in zip(expected_annotations, parsed_annotations):
        assert int(expected['start']) == parsed['start']
        assert int(expected['end']) == parsed['end']
        assert int(expected['strand']) == parsed['strand']
        assert expected['category'] == parsed['dnafeature']['category']
        assert expected['name'] == parsed['dnafeature']['name']
        assert int(expected['length']) == parsed['dnafeature']['length']

def test_parse_genbank_unrecognized_feature_type():
    """Can parse a genbank file with an unrecognized feature; feature is retained"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/unrecognized_feature_type.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            print data
            assert data['locus']
            assert len(data['dnafeatures']) == 1


def test_parse_genbank_unrecognized_feature_qualifier():
    """Can parse a genbank file with an unrecognized feature qualifier; qualifier is retained"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/unrecognized_feature_qualifier.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            print data
            assert data['locus']
            assert len(data['dnafeatures']) == 1
            annotation = data['dnafeatures'][0]
            assert annotation['dnafeature']['properties']['unrecognizable_feature_qualifier']


def test_parse_genbank_unrecognised_missing_feature_header():
    """Can parse a genbank file with missing feature header"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/missing_feature_header.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        parsed_annotations = ret['dnafeatures']


def test_parse_genbank_no_features():
    """Can parse a genbank file with feature header but no features"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/no_features.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            assert data['locus']
            assert len(data['dnafeatures']) == 0


def test_parse_genbank_missing_one_feature():
    """Can parse a genbank file with one feature"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/one_feature.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            assert data['locus']
            assert len(data['dnafeatures']) == 1
            annotation = data['dnafeatures'][0]
            assert annotation['dnafeature']['pattern']['bases'] == 'CCCCTT'


def test_parse_genbank_ambiguous_dna():
    """Can parse a genbank file with ambiguous DNA chars in the pattern"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/ambiguous.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            assert data['locus']
            assert len(data['dnafeatures']) == 1
            annotation = data['dnafeatures'][0]
            assert annotation['dnafeature']['pattern']['bases'] == 'CCMCTT'


def test_parse_single_base_coordinate():
    """Can parse feature coordinates which are only one base long"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/single_base_feature.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            assert data['locus']
            assert len(data['dnafeatures']) == 1
            annotation = data['dnafeatures'][0]
            assert annotation['dnafeature']['length'] == 1


def test_feature_spans_plasmid_coordinate_origin():
    """Can parse feature spanning the origin of plasmid circular coordinate system"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../../data/genbank/span_circular.gb')
    with open(test_data_path, 'rb') as input_fh:
        ret = genbank.parse(input_fh)
        for data in ret if isinstance(ret, list) else [ret]:
            assert data['locus']
            assert len(data['dnafeatures']) == 1
            annotation = data['dnafeatures'][0]
            assert annotation['dnafeature']['length'] == 20
            assert int(annotation['start']) == 20
            assert int(annotation['end']) == 10
            assert annotation['dnafeature']['pattern']['bases'] == 'GGGGGGGGGGAAAAAAAAAA'
