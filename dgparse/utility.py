#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Utility methods for working with sequence parsers
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals


import re
from string import maketrans

# Anything that is not DNA
NOT_AMBIG_DNA = re.compile(r"[^ACGTacgtMmRrWwSsYyKkVvHhDdBbXxNn]")

NOT_DNA = re.compile(r"[^ACGTacgt]") #anything that's not ACGT

DNA_CHAR = set("ACGT")

trans_table = maketrans(b"ACGTacgtMRWSYKVHDBXNmrwsykvhdbxn",
                        b"TGCATgcaKYWSRMBDHVXNkywsrmbdhvxn")


def get_complement(seq_str):
    '''get the complement of this sequence G --> C'''
    #Sequence will always be DNA in this context
    normalised_ascii = seq_str.encode('ascii', 'ignore')


    return normalised_ascii.translate(trans_table)


def get_reverse_complement(seq_str):
    '''complement a sequence and reverse it to get the rc sequence that binds
    to the input sequence'''
    compliment = get_complement(seq_str)
    return compliment[::-1]