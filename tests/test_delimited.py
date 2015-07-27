"""
Unit tests for the delimited file (CSV) parser.
"""

import os
import pytest

from dgparse import delimted

@pytest.fixture(params=[
    '../data/delimited/features.csv',
    '../data/delimited/sequences.csv'
])
def record_buffer(request):
    """Contains an OPEN test file or buffer"""
    test_file_path = os.path.join(os.path.dirname(__file__), request.param)
    test_file = open(test_file_path, 'r')

    def teardown():
        test_file.close()

    request.addfinalizer(teardown)
    return test_file


def test_csv_parse(record_buffer):
    """Test invokation of the parser"""
    records = []
    records.extend(delimted.parse(record_buffer))  # path or buffer?
    assert len(records) > 0


