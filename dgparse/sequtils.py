# encoding=utf-8
"""
A collection of functions for dealing with sequences. Generally designed
to accept an input string and out either a modified output string or list of
results.
"""
from __future__ import unicode_literals, division, absolute_import

import hashlib
import re
import string

TRANSLATE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
    'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
    'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L',
    'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I',
    'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T',
    'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A',
    'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D',
    'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G', 'TAA': "*", 'TAG': "*", 'TGA': "*",
}

AA_SHORTNAME = {
    'G': 'Gly',
    'P': 'Pro',
    'A': 'Ala',
    'V': 'Val',
    'L': 'Leu',
    'I': 'Ile',
    'M': 'Met',
    'C': 'Cys',
    'F': 'Phe',
    'Y': 'Tyr',
    'W': 'Trp',
    'H': 'His',
    'K': 'Lys',
    'R': 'Arg',
    'Q': 'Gln',
    'N': 'Asn',
    'E': 'Glu',
    'D': 'Asp',
    'S': 'Ser',
    'T': 'Thr',
    '*': '*',
}


DNA_COMPLEMENTS = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
}

ASCII_DNA_COMP = {str(k): str(v) for k, v in DNA_COMPLEMENTS.iteritems()}

AMBIGUOUS_CODES = {
    'W': ['A', 'T'],
    'N': ['A', 'C', 'G', 'T'],
    'S': ['C', 'G'],
    'M': ['A', 'C'],
    'K': ['G', 'T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
}


DNA_BYTES = {"A": b"0", "C": b"1", "G": b"2", "T": b"3"}

DNA_CHAR = set("ACGT")
AMBIG_CHAR = set("ACGTacgtMmRrWwSsYyKkVvHhDdBbXxNn")
NOT_UNAMBIG_DNA = re.compile(r"[^ACGTacgt]") #anything that's not ACGT
MOD_CHAR = set("*")

#Anything that is not
NOT_DNA = re.compile(r"[^ACGTacgtMmRrWwSsYyKkVvHhDdBbXxNn]")

UNICODE_TABLE = dict((ord(key), value) for key, value in
                DNA_COMPLEMENTS.iteritems())


def get_complement(seq_str):
    """get the complement of this sequence G --> C"""
    #Sequence will always be DNA in this context
    if isinstance(seq_str, str):
        trans_table = string.maketrans(b'ACGTacgt', b'TGCAtgca')
        return seq_str.translate(trans_table)
    else:
        return seq_str.translate(UNICODE_TABLE)


def get_reverse(seq_str):
    """reverse the sequence, read 3' to 5'  """
    output = seq_str[::-1]
    return output


def get_reverse_complement(seq_str):
    """complement a sequence and reverse it to get the rc sequence that binds
    to the input sequence"""

    compliment = get_complement(seq_str)
    return compliment[::-1]


def compute_sha1(data):
    """Compute the sha1 hash of a sequence"""
    try:
        bases = data.get('bases').replace('\n', '').upper()
        data['sha1'] = hashlib.sha1(bases).hexdigest()
        data['bases'] = bases
    except AttributeError:
        pass  # let the validator catch this
    finally:
        return data


rev_comp = get_reverse_complement
comp = get_complement


def dotsetter(tokens, value, result):
    """Set an attribute in a nested JSON blob"""
    key = tokens.pop()
    if tokens:
        if key not in result:
            result[key] = dict()
        dotsetter(tokens, value, result[key])
    else:
        result.update({key: value})
