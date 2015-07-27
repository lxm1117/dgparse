# -*- coding: utf-8 -*-
"""
Unit tests for the delimited file (CSV) parser.
"""

import os
import pytest

from dgparse import delimited
from dgparse import exc

@pytest.fixture(params=[
    '../data/delimited/dnafeature.csv',
    '../data/delimited/plasmid.csv'
])
def record_buffer(request):
    """Contains an OPEN test file or buffer"""
    test_file_path = os.path.join(os.path.dirname(__file__), request.param)
    test_file = open(test_file_path, 'r')

    def teardown():
        test_file.close()

    request.addfinalizer(teardown)
    return test_file


@pytest.fixture
def nonrecordline_csv(request):
    """
    A CSV file with
    :return:
    """
    path = '../data/delimited/non_record_line.csv'
    test_file_path = os.path.join(os.path.dirname(__file__), path)
    test_file = open(test_file_path, 'r')

    def teardown():
        test_file.close()

    request.addfinalizer(teardown)
    return test_file


def test_csv_parse(record_buffer):
    """Test invokation of the parser"""
    records = []
    records.extend(delimited.parse(record_buffer))  # path or buffer?
    for record in records:
        if record['type_'].lower() == 'dnafeature':
            assert 'pattern' in record
        elif record['type_'].lower() == 'plasmid':
            assert 'sequence' in record


