#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Module contains functions for looking up or guessing the file format of
various sequence files.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import genbank
import snapgene
import fasta

EXTENSIONS = {
    '.gb': genbank.parse,
    '.gbk': genbank.parse,
    '.genbank': genbank.parse,
    '.fa': fasta.parse_fasta_file,
    '.fasta': fasta.parse_fasta_file,
    '.dna': snapgene.main.parse,
}


def get_parser(extension):
    """get or guess the parser for an file"""
    try:
        return EXTENSIONS[extension]
    except KeyError:
        raise

