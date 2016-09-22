#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Unit tests for schema adaptation and validation behaviour
"""

import pytest
import json
from dgparse import schema
import uuid

@pytest.mark.parametrize("name,bases_in,bases_out,expected_errs,expected_mods", [
    ('pass',     'AGTCAGTCAGTC',   'AGTCAGTCAGTC',   {}, []),
    ('space',    'AGTCAGT CAGTC',  'AGTCAGT CAGTC',  {'sequence': {'bases': ['Non-IUPAC DNA base found at 7']}}, []),
    ('dash',     'AGTCAGT-CAGTC',  'AGTCAGT-CAGTC',  {'sequence': {'bases': ['Non-IUPAC DNA base found at 7']}}, []),
    ('ambig',    'AGTCAGTNCAGTC',  'AGTCAGTNCAGTC',  {}, []),
    ('modone',   'AGTCAGT*CAGTC',  'AGTCAGTCAGTC',   {}, [{'position': 7, 'symbol': '*'}]),
    ('modtwo',   'AGTCAGT*C*GTCC', 'AGTCAGTCGTCC',   {}, [{'position': 7, 'symbol': '*'}, {'position': 9, 'symbol': '*'}]),
    ('dashmod',  'AGTCAGT*C-GTC',  'AGTCAGTC-GTC',   {'sequence': {'bases': ['Non-IUPAC DNA base found at 8']}}, [{'position': 7, 'symbol': '*'}]), 
])
def test_sequence_schema(name, bases_in, bases_out, expected_errs, expected_mods):
    '''Parser identifies and cleans non-standard chars as expected'''
    raw = {
        'accession': uuid.uuid4().hex,
        'name': name,
        'sequence': {'bases': bases_in}
    }
    validator = schema.DnaOligoSchema()
    loaded, errors = validator.load(raw)
    json.dumps(loaded, indent=4)
    assert errors == expected_errs
    assert loaded['modifications'] == expected_mods
    if not expected_errs:
        # In marshmallow 2.9.1 erroneous fields are NOT returned
        # This should have always been the case.
        # Therefore the original data not the validated object should be inspected
        # for errors.
        assert loaded['sequence']['bases'] == bases_out
