#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Centralized Parser for models coming from Microsoft Excel Spreadsheets
"""
from functools import partial
import openpyxl


def row_to_dict(headers, row_data):
    """Convert a row to a valid dictionary"""
    record = {key: value for key, value in zip(headers, row_data)}
    return record


def make_excel_parser(type_):
    """Constructor which returns an excel parser callable for a given model"""

    def excel_parser(xls_file_path):
        """Takes a file object and returns a list of records"""
        records = []
        with openpyxl.load_workbook(xls_file_path, read_only=True) as wbook:
            wsheet = wbook[type_]  # name of sheet is always the record type
            headers = wsheet.rows[0]  # always the attributes
            record_factory = partial(row_to_dict, headers)
            records.extend(map(record_factory(wsheet.rows[1:])))
        return records

    return excel_parser