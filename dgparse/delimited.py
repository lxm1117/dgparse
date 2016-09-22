#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Centralized Parser for models coming from delimiter separated value files
"""

import csv
import os
import functools

from .sequtils import dotsetter


def clean_record(basename, record):
    """Clean the extracted record and nest fields"""
    result = {'__class__': basename.split('.')[0]}
    for key, value in record.iteritems():
        value = '' if value is None else value
        value = True if value in {'True', 'true', 'TRUE'} else value
        value = False if value in {'False', 'false', 'FALSE'} else value
        if key is '':
            msg = "{0} contains a NonRecord Entry {1}".format(basename, record)
            result['ERROR'] = msg
        if '.' in key:
            tokens = key.split('.')
            tokens.reverse()  # we pop tokens so this is critical
            dotsetter(tokens, value, result)
        else:
            result[key] = value
    return result


def parse(open_file, fieldnames=None, record_type=None, delimiter=b","):
    """
    Parse an open file object
    :param open_file:
    :return:
    """
    # What does this file contain?
    if not record_type:
        record_type = os.path.basename(getattr(open_file, 'name', 'na.unk'))
    record_cleaner = functools.partial(clean_record, record_type)
    return [record_cleaner(row)
            for row in csv.DictReader(open_file, fieldnames,
                                      delimiter=delimiter)]
