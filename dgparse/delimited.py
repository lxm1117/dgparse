#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Centralized Parser for models coming from delimiter separated value files
"""

import csv
import os

from dgparse import exc

def clean_record(record):
    if '' in record:  # blank key
        msg = "{0} contains a NonRecord Entry {1}".format(open_file.name, record)
        raise exc.NonRecordEntry(msg)
    record.update({"type_": record_type})
    else:
        return record

def parse(open_file, fieldnames=None, record_type=None):
    """
    Parse an open file object
    :param open_file:
    :return:
    """
    records = []  # always return a list
    # What does this file contain?
    if not record_type:
        record_type = os.path.basename(open_file.name)

    for record in csv.DictReader(open_file, fieldnames):


        records.append(record)
    return records
