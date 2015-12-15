# -*- coding: utf-8 -*-
"""
Unit tests for the delimited file (CSV) parser.
"""

import os
import io
import pytest
import collections

from dgparse import fasta
from dgparse import exc

@pytest.mark.parametrize("path,expected_name,expected_length", (
    ('upper_case.fasta', 'upper_case', 16),
    ('upper_case_multiline.fasta', 'upper_case_multiline', 16),
))
def test_parse_fasta(path, expected_name, expected_length):
    """Can parse file and retrieve expected name and length"""
    test_file_path = os.path.join(os.path.dirname(__file__), '../../data/fasta/' + path)
    with open(test_file_path, 'rb') as test_file:
        ret = fasta.parse(test_file)
        ret_dict = ret if isinstance(ret, dict) else ret[0] # first sequence in file
        assert ret_dict['name'] == expected_name
        assert len(ret_dict['sequence']['bases']) == expected_length
