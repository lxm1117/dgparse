#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Unit tests for schema adaptation and validation behaviour
"""

import pytest
import uuid
import copy
from dgparse import schema

TEST_BASES = 'ATCGATCGATCGATCGATCGATCGATCG'
TEST_SEQUENCE = {
    'accession': uuid.uuid4().hex,
    'name': 'test',
    'sequence': {'bases': TEST_BASES},
}
TEST_PATTERN = {
    'accession': uuid.uuid4().hex,
    'name': 'test',
    'pattern': {'bases': TEST_BASES},
}

@pytest.mark.parametrize("length", [
    None,
    len(TEST_BASES),
    len(TEST_BASES) + 1,
    0,
])
def test_length_with_specified_length(length):
    """Length is calculated"""
    for input_data, validator in [
            (TEST_SEQUENCE, schema.DnaOligoSchema()),
            (TEST_SEQUENCE, schema.DnaMoleculeSchema()),
            (TEST_PATTERN, schema.DnaFeatureSchema()),
    ]:
        local_data = copy.deepcopy(input_data)
        local_data.update({'length': length})
        loaded, errors = validator.load(local_data)
        assert loaded['length'] == len(TEST_BASES)
        assert errors == {}

def test_length_without_specified_length():
    """Length is calculated"""
    for input_data, validator in [
            (TEST_SEQUENCE, schema.DnaOligoSchema()),
            (TEST_SEQUENCE, schema.DnaMoleculeSchema()),
            (TEST_PATTERN, schema.DnaFeatureSchema()),
    ]:
        loaded, errors = validator.load(input_data)
        assert loaded['length'] == len(TEST_BASES)
        assert errors == {}
