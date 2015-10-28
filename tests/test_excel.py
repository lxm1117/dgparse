"""
Unit tests for the excel parser.
"""

import os
import pytest
import xlsxwriter

from dgparse import excel


@pytest.fixture(params=[
    '../data/excel/oligos.xlsx',
    '../data/excel/primers.xlsx'
])
def record_buffer(request):
    """Contains an OPEN test file or buffer"""
    test_file_path = os.path.join(os.path.dirname(__file__), request.param)
    test_file = open(test_file_path, 'rb')

    def teardown():
        test_file.close()

    request.addfinalizer(teardown)
    return test_file


def test_excel_parse(record_buffer):
    """Test invokation of the parser"""
    records = []
    records.extend(excel.parse(record_buffer))  # path or buffer?
    assert len(records) > 0


def test_compound_keys_removed(record_buffer):
    """Compound keys containing a '.' e.g. sequence.bases are removed"""
    records = []
    records.extend(excel.parse(record_buffer))
    for record in records:
        for key in record.keys():
            assert '.' not in key

