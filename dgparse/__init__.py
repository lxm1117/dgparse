#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Desktop Genetics Parsing and Valdiation Python Package
This package contains parsers for a wide range of biology file formats
and accompanying schema validators and serialization code to ensure
compatibility with the current DeskGen genome editing platform domain model.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import os
import logging
import functools

from . import schema
from . import exc
from . import delimited
from . import snapgene
from . import excel
from . import genbank
from . import fasta

VALIDATORS = {
    'oligo': schema.DnaOligoSchema(),
    'primer': schema.DnaPrimerSchema(),
    'plasmid': schema.DnaPlasmidSchema(),
    'construct': schema.DnaConstructSchema(),
    'dnafeature': schema.DnaFeatureSchema(),
    'dnamolecule': schema.DnaMoleculeSchema(),
    'dnadesign': schema.DnaDesignSchema(),
}

PARSERS = {
    '.csv': delimited.parse,
    '.tsv': functools.partial(delimited.parse, delimiter=b"\t"),
    '.dna': snapgene.parse,
    '.xlsx': excel.parse,
    '.gb': genbank.parse,
    '.gbk': genbank.parse,
    '.genbank': genbank.parse,
    '.fa': fasta.parse,
    '.fasta': fasta.parse,
}

LOG = logging.getLogger(__file__)


def validate(record):
    """
    Returns
    :param type_:
    :param record:
    :return:
    """
    if 'ERROR' in record:
        raise exc.FormatException(record['ERROR'])
    type_ = record.pop('__class__')
    try:
        validator = VALIDATORS[type_]
    except KeyError:
        msg = "No record type defined for {0}".format(record)
        raise exc.UndefinedRecordType(msg)
    # IMPORTANT: load doesn't construct an object but MAPS it to the new schema
    data, errors = validator.load(record)
    return data, errors


def load_iter(record_type, record_files):
    """Load a file and return a single array of records"""
    record_schema = VALIDATORS[record_type]
    for record_path in record_files:
        format_ = os.path.splitext(record_path)[-1]
        parser = PARSERS[format_]
        with open(record_path, 'r') as record_file:
            try:
                raw_records = parser(record_file)
                if isinstance(raw_records, dict):
                    yield record_schema.load(raw_records)
                else:
                    for item in raw_records:
                        yield record_schema.load(item)
            except exc.ParserException as exception:
                LOG.error(exception)
                raise exception
