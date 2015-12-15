# -*- coding: utf-8 -*-
"""
Test genbank feature annotations
"""

import os
import io
import csv
import pytest

from dgparse import genbank

def test_parse_fasta():
    """Can parse genbank file and retrieve expected features"""

    test_data_path = os.path.join(os.path.dirname(__file__), '../data/genbank/PX330.gbk')
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
