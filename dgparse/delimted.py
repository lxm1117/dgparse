#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Centralized Parser for models coming from delimiter separated value files
"""

import csv


def parse(open_file, fieldnames=None):
    """
    Parse an open file object
    :param open_file:
    :return:
    """
    # stupid but maintains a uniform interface
    return csv.DictReader(open_file, fieldnames)
