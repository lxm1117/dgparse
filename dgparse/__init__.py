#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script does something ...
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

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
    'dnafeature': schema.DnaFeatureSchema(),
    'dnamolecule': schema.DnaMoleculeSchema(),
    'dnadesign': schema.DnaDesignSchema(),
}

PARSERS = {
    '.csv': delimited.parse,
    '.dna': snapgene.parse,
    '.xlsx': excel.parse,
    '.gb': genbank.parse,
    '.gbk': genbank.parse,
    '.genbank': genbank.parse,
    '.fa': fasta.parse,
    '.fasta': fasta.parse,
}


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
