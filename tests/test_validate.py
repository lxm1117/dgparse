"""
Test validation that parsed data matches DeskGen Schema
"""
from __future__ import unicode_literals, division
import os
import dgparse


def parse_record_array(file_path):
    """Returns a PythonList of PyDict records"""
    test_file_path = os.path.join(os.path.dirname(__file__), file_path)
    path, extension = os.path.splitext(test_file_path)
    parser = dgparse.PARSERS[extension]
    records = []
    with open(test_file_path, 'r') as record_buffer:
        records.extend(parser(record_buffer))

    return records


def test_validate_plasmids():
    record_array = parse_record_array('../data/delimited/plasmid.csv')
    for number, record in enumerate(record_array):
        data, errors = dgparse.validate(record)
        assert data['length'] > 0
        assert errors == {}
        assert len(data['sequence']['bases']) == int(data['length'])


def test_weird_char():
    record_array = parse_record_array('../data/delimited/dnafeature.csv')
    weird_record = record_array[7]
    data, errors = dgparse.validate(weird_record)
    assert errors


def test_validate_features():
    record_array = parse_record_array('../data/delimited/dnafeature.csv')
    for number, record in enumerate(record_array):
        data, errors = dgparse.validate(record)
        if errors is {}:
            assert data['length'] > 0
            assert len(data['pattern']['bases']) == data['length']
