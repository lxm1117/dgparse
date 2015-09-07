# encoding=utf-8
"""
Genbank Schema Constants
"""
GENBANK_HEADERS = (
    # standard
    'LOCUS',
    'DEFINITION'
    'ACCESSION'
    'VERSION'
    'KEYWORDS'
    'SOURCE'
        #'ORGANISM'
    'REFERENCE'
        #'AUTHORS'
        #'TITLE'
        #'JOURNAL'
        #'PUBMED'
    'FEATURES',
    'ORIGIN',

    # found
    'BASE COUNT',
)
GENBANK_DIVISONS = (
    'PRI',  # primate sequences
    'ROD',  # rodent sequences
    'MAM',  # other mammalian sequences
    'VRT',  # other vertebrate sequences
    'INV',  # invertebrate sequences
    'PLN',  # plant, fungal, and algal sequences
    'BCT',  # bacterial sequences
    'VRL',  # viral sequences
    'PHG',  # bacteriophage sequences
    'SYN',  # synthetic sequences
    'UNA',  # unannotated sequences
    'EST',  # expressed sequence tags
    'PAT',  # patent sequences
    'STS',  # sequence tagged sites
    'GSS',  # genome survey sequences
    'HTG',  # high-throughput genomic sequences
    'HTC',  # unfinished high-throughput cDNA sequencing
    'ENV',  # environmental sampling sequences
    'UNK',  # unknown
)
VALID_ORIENTATIONS = ('linear', 'circular')
ORIENTATION_MAPPING = dict(zip(VALID_ORIENTATIONS, [False, True]))
KNOWN_DATE_FORMATS = ('%d-%b-%Y', '%d-%m-%Y')
