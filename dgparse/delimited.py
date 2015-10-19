#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Centralized Parser for models coming from delimiter separated value files
"""

import csv
import os
import functools

def clean_record(basename, record):
    result = {'__class__': basename.split('.')[0]}
    for key, value in record.iteritems():
        value = '' if value is None else value.replace(' ', '')
        if key is '':
            msg = "{0} contains a NonRecord Entry {1}".format(basename, record)
            result['ERROR'] = msg
        if '.' in key:
            parent_attr, child_attr = key.split('.')
            if parent_attr in result:
                result[parent_attr][child_attr] = value
            else:
                result[parent_attr] = {child_attr: value}
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
