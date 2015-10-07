# -*- coding: utf-8 -*-
"""
Unit tests for the delimited file (CSV) parser.
"""

import os
import io
import pytest
import collections

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
        if record['__class__'].lower() == 'dnafeature':
            assert 'pattern' in record
        elif record['__class__'].lower() == 'plasmid':
            assert 'sequence' in record


@pytest.fixture()
def upload_request():
    """A dummy HTTP upload request"""
    file_dict = {
        'file': io.StringIO(
            u'name,sequence.bases\nseq1,AGTCAGTCAGTCAGTCAGTC\n'),
        'filename': 'thing',
    }
    file_ = collections.namedtuple(
        'FileDict', file_dict.keys())(**file_dict)
    upload_dict = {
        'POST': {'upload': file_},
    }
    upload_ = collections.namedtuple(
        'UploadDict', upload_dict.keys())(**upload_dict)
    return upload_


def test_csv_parse_http_upload(upload_request):
    """
    Parse an HTTP upload request
    """
    ret = delimited.parse(upload_request.POST.get('upload').file)
    for record in [ret] if isinstance(ret, dict) else ret:
        assert 'name' in record
        assert record['name'] == 'seq1'
