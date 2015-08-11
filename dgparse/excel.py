#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Centralized Parser for models coming from Microsoft Excel Spreadsheets
"""

from functools import partial
import logging

import openpyxl

log = logging.getLogger(__file__)


def row_to_dict(headers, constants, row_data):
    """Convert a row to a valid dictionary"""
    record = {}
    record.update(constants)
    if row_data is None:  # handle dead rows
        return
    for key, value in zip(headers, row_data):
        if key.value is None:
            continue
        if hasattr(key, 'value'):
            key = key.value
        if hasattr(value, 'value'):
            value = value.value
        if '.' in key:
            key, child_attr = key.split('.')
            value = {child_attr: value}
        record[key] = value
    return record


def parse(open_file):
    """Constructor which returns an excel parser callable for a given model"""
    records = []  # Use empty list to represent an empty set
    wbook = openpyxl.load_workbook(open_file)  #turn zipped file into workbook
    sheets = wbook.get_sheet_names()  # get the types of records included
    for sheet in sheets:
        wsheet = wbook[sheet]  # name of sheet is always the record type
        headers = wsheet.rows[0]  # always the attributes
        record_factory = partial(row_to_dict, headers, {})  # no constants yet
        records.extend(filter(None, map(record_factory, wsheet.rows[1:])))

    return records