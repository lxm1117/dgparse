#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Unit tests for schema adaptation and validation behaviour
"""

from __future__ import unicode_literals, division, absolute_import
import os
import pytest
from dgparse import excel
from dgparse import schema


@pytest.fixture(scope='session')
def oligo_collection():
    oligo_path = os.path.join(os.path.dirname(__file__), '../data/excel/oligos.xlsx')
    with open(oligo_path, 'r') as oligo_file:
        raw_records = excel.parse(oligo_file)
    # the second column is the test case name
    oligos = {r['name']: r for r in raw_records}
    return oligos


@pytest.fixture(scope='session')
def primer_collection():
    primer_path = os.path.join(os.path.dirname(__file__), '../data/excel/primers.xlsx')
    with open(primer_path, 'r') as primer_file:
        raw_records = excel.parse(primer_file)
    primers = {r['name']: r for r in raw_records}
    return primers


@pytest.mark.parametrize("case_name,expected", [
    ('good_oligo', 'accept'),
    ('mixed_case', 'accept'),
    ('modified_oligo', 'accept'),
    ('no_delta_g', 'accept'),
    ('lower_case', 'accept'),
    ('no_target', 'accept'),
    ('lower_mods', 'accept'),
    ('mixed_mods', 'accept'),
    ('no_data', 'reject'),
    ('no_accession', 'reject')
])
def test_load_oligo(case_name, expected, oligo_collection):
    """Test validation of oligos"""
    oligo_schema = schema.DnaOligoSchema()
    input_data = oligo_collection[case_name]
    oligo, errors = oligo_schema.load(input_data)
    if expected == 'accept':
        assert oligo
        assert errors == {}
    else:
        assert errors


@pytest.mark.parametrize("case_name,expected", [
    ('no_data_f', 'reject'),
    ('no_data_r', 'reject'),
    ('x_conc_f', 'accept'),
    ('x_conc_r', 'accept'),
    ('split_units', 'accept'),
    ('priming_seq', 'accept'),
    ('homology_seq', 'accept'),
    ('lower_seq', 'accept'),
    ('good', 'accept'),
    ('mixed_case', 'accept')
])
def test_load_primer(case_name, expected, primer_collection):
    """Test validation of primers"""
    primer_schema = schema.DnaPrimerSchema()
    input_data = primer_collection[case_name]
    primer, errors = primer_schema.load(input_data)
    if expected == 'accept':
        assert primer
        assert errors == {}
    else:
        assert errors
