Desktop Genetics Parsing Library
================================

So many file formats, so little time

## Purpose and Design
The dgparse Python library contains several utility functions for transforming
one sort of DNA file format into another. While other libraries strive for 
strict enforcement of a particular standard, this library's parsers aim to be as
resilient as possible given the enormous variety of file formats in use today.
For example, a basic parser may extract only the sequence from a file if the 
metadata is corrupt or indecipherable, but a more advanced parser would also 
extract some of this metadata. Parsers should be *greedy* and try to extract as
much useful data as possible from a file, even if it is incorrectly formatted.


## Interface
All parsers in this library **must** implement a common interface. The Input
to the parser shall be a "buffer-like" (file-like) Object. A utility 
function or mapping will have already attempted to determine the file type by 
it's extension or other clues. An optional Input shall be a dictionary of default
values.

Example:

```
with open('target_file_path', 'r') as target_file_input:
    output = parse_file(target_file_input, defaults={})
```

The output of the parser shall be a JSON-serializable set of Python Primitives
(Lists, Dictionaries, Strings, Integers, etc). This format largely parallels that
of BioPython, except it is readily serialized to JSON or YAML and thus very 
portable between systems, libraries, and languages.

Example:

```
output = [ # This outer array accomidates multi-record fasta files and similar formats
    {
        'name': 'pBR322',
        'accession': 'A universially unique identifier',
        'namespace': 'Where that unique identifier comes from',
        'description': 'Text Blob',
        'is_circular': True or False,
        'properties': {},  # assorted key value pairs
        'sequence': {
            'bases': "ACGT",
            'sha1': "The 40 char sha1.hexdigest of the sequence",
        },
        'dnafeatures': [
        # array of dnafeatures
            {
                'name': "Name of the DNA Feature",
                'start': 0, # in the sequence
                'end': 1, # integer start and end positions,
                'strand': -1, 0, 1 # -1 is the compliment strand
                'feat_type': 'the type of this feature',
                'pattern': 'Bases of the DNA',
                'description': 'description of the feature',
                'properties': 'as before, but for this feature'
            },
        ],
        'restrictionsites': [
        # array of restriction sites, the same as features
        ]
    }
# other records
]
```
## Exporting
To export to a specific file format the output object described above is used
as the input to a template file. Python is full of rich templating languages and
we have found this to be a very extensible way of outputting an exotic format.
Otherwise, one can simply use the standard JSON, YAML, or msgpack libraries to
obtain very compact, portable serializations of the above object.