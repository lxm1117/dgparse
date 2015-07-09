
"""
A collection of functions for dealing with sequences. Generally designed
to accept an input string and out either a modified output string or list of
results.

Imported from AutoClone
"""
from __future__ import unicode_literals

import itertools
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
    'GGG': 'G', 'TAA': "*", 'TAG': "*", 'TGA': "*", }

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

ambiguous_codes = {'W': ['A', 'T'],
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

NOT_UNAMBIG_DNA = re.compile(r"[^ACGTacgt]") #anything that's not ACGT

#Anything that is not
NOT_DNA = re.compile(r"[^ACGTacgtMmRrWwSsYyKkVvHhDdBbXxNn]")

unicode_table = dict((ord(key), value) for key, value in
                DNA_COMPLEMENTS.iteritems())

GC_CONTENT_REGEX = re.compile(r'[GC]')
CONSECUTIVE_BASE_REGEX = re.compile(r'(.)\1{3,}')
SAPI_FWD = 'GCTCTTC'
SAPI_REV = 'GAAGAGC'


def get_complement(seq_str):
    """get the complement of this sequence G --> C"""
    #Sequence will always be DNA in this context
    if isinstance(seq_str, str):
        trans_table = string.maketrans(b'ACGTacgt', b'TGCAtgca')
        return seq_str.translate(trans_table)
    else:
        return seq_str.translate(unicode_table)


def reverse(seq_str):
    """reverse the sequence, read 3' to 5'  """
    output = seq_str[::-1]
    return output


def get_reverse_complement(seq_str):
    """complement a sequence and reverse it to get the rc sequence that binds
    to the input sequence"""

    compliment = get_complement(seq_str)
    return compliment[::-1]

rev_comp = get_reverse_complement
comp = get_complement

def compare_sequences(seq1, seq2):
    """Ideally an implementation of the SW or NW alignment algorithms"""
    #print len(seq1), len(seq2)
    for i, base in enumerate(seq1):
        if i > len(seq2):
            break
        if i<len(seq2):
            print i, base, seq2[i]
            if seq2[i] != base:
                print "FAIL at", i, base, seq2[i]


class Keeper(object):
    def __init__(self, keep):
        self.keep = set(map(ord, keep))

    def __getitem__(self, n):
        if n not in self.keep:
            return None
        return unichr(n)

    def __call__(self, s):
        return unicode(s).translate(self)


makefilter = Keeper


def translate(seq, table):
    seq = seq.upper()
    n = len(seq)
    amino_acids = [table[seq[i:i+3]]for i in range(0, n-n%3, 3)]
    return "".join(amino_acids)


def find_mismatches(left, right):
    mismatches = []
    for i, c in enumerate(left):
        if right[i] != c:
            mismatches.append(i)
    return mismatches


class SequenceCoordinates(object):
    """Generalized Sequence coordinates for use outside DB env"""
    __slots__ = ('chromosome_id', 'start', 'end', 'strand')

    def __init__(self, chromosome_id, start_end, strand):
        # we split inside but keep the compatible interface
        self.chromosome_id = chromosome_id
        try:
            self.start = start_end[0]  # a Tuple or Numeric Range
            self.end = start_end[1]
        except (TypeError, AttributeError):  # The error raised is unfortunately not consistent
            self.start = start_end.lower
            self.end = start_end.upper
        self.strand = strand

    def __eq__(self, other):
        """
        To coordinates are equal if the specify the same place on the same
        chromosome of the same genome
        :param other:
        :return:
        """
        return all((
            (self.chromosome_id == other.chromosome_id),
            (self.start == other.start),
            (self.end == other.end),
            (self.strand == other.strand))
        )

    def __hash__(self):
        """
        Make coordinates hashable so they can be compared to each other.
        Necessary for asserting if two things have the same coordinates.
        :return:
        """
        # We can split out the Numeric Range Type if needed
        # Chromosome id is used as we must be sure we're talking about a
        # precise sequence, not just any copy of chrY or chrX
        return hash((self.chromosome_id, self.start, self.end, self.strand))

    def __iter__(self):
        return (self.chromosome_id,
                self.start,
                self.end,
                self.strand).__iter__()

    def __repr__(self):
        return "{0.chromosome_id} {0.start} {0.end} {0.strand}".format(self)

    def serializable(self):
        return {
            'chromosome_id': self.chromosome_id,
            'start_end': (int(self.start), int(self.end)), # don't break the interface
            'strand': self.strand
        }

    def __json__(self, request):
        return self.serializable()

    @property
    def start_end(self):
        """Support the standard coordinates object more cleanly"""
        return self.start, self.end


def bases_codons(bases, phase=0):
    'Return codons in a sequence according to phase'
    # uses ensembl "standard" to interpret what you mean by phase
    # http://lists.ensembl.org/ensembl-dev/msg00562.html
    # http://upload.wikimedia.org/wikipedia/en/d/db/Exon_and_Intron_classes.png
    phase_slice = (3 - phase) % 3
    it_bases = iter(bases[phase_slice:])
    while True:
        chunk = ''.join(itertools.islice(it_bases, 3))
        if not chunk or len(chunk) != 3:
            return
        yield chunk


def codons_polypeptides(codons):
    'Return a string of polypeptides for a list of codons'
    return ''.join([TRANSLATE[codon] for codon in codons])


def translate_sequence(bases):
    'Return alternative translations for a sequence of basepairs by phase'
    return [
        {
            'phase': phase,
            'polypeptides': codons_polypeptides(bases_codons(bases, phase))
        }
        for phase in range(3)
    ]


def case_diff(bases, alternative):
    'Return alternative sequence with differing bases in lowercase'
    for index, base in enumerate(alternative):
        if base == bases[index]:
            yield base
        else:
            yield base.lower()
