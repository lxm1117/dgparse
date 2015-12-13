# -*- coding: utf-8 -*-
"""
Unit tests for the genbank parser.
"""

import os
import io
import csv
import pytest
import collections
import json

from dgparse import genbank
from dgparse import exc

@pytest.mark.parametrize("path,", (
    '../data/genbank/PX330.gbk',
))
def test_parse_fasta(path):
    """Can parse file and retrieve expected features"""
    test_file_path = os.path.join(os.path.dirname(__file__), path)
    with open(test_file_path, 'rb') as test_file:
        ret = genbank.parse(test_file)
        parsed_names = [o['dnafeature']['name'] for o in ret['dnafeatures']]
        expected_names = [
            'CBh',
            'hSpCsn1',
            'bGH polyA',
            'R-ITR',
            'f1 Origin',
            'AmpR',
            'pUC ori',
            'Human U6 Promoter',
            'NLS',
            'chimeric guide RNA scaffold',
            'U6 terminator',
            '3xFLAG',
            'NLS',
            'EcoRI',
            'T2A',
            'GFP',
            'EcoRI'
        ]
        assert len(parsed_names) == len(expected_names)
        assert parsed_names == expected_names
