"""Parse GTF for Gene, Transcript, Exon and CDS features

Usage:
  parse_gtf.py (--input=<file.gz>) (--template=<file>) (--config=<file>) (--namespace=<str>) [--output=<dir>] [--require-ccds=<bool>] [--translate=<bool>] [--transform=<bool>] [--log=<loglevel>]

Options:
  --help                    Show this screen
  --input=<file>            GTF file, gzipped
  --template=<file>         JSON template for I/O field mapping
  --config=<file>           Genome config, CSV file specifying the seqnames to extract features for and other data
  --namespace=<str>         Default name space
  --transform=<bool>        Prepare files ready for database loading [default: 1]
  --require-ccds=<bool>     Output cds/cdsregions only if there is a CCDS identifier [default: 1]
  --translate=<bool>        Attempt to construct and translate the coding sequence [default: 0]
  --log=<loglevel>          Logging level [default: ERROR]
  --output=<dir>            Output/working directory [default: ./]

"""
import os
import re
import csv
import gzip
import logging
import collections
import json
import yaml

#from bio.sequtils import get_reverse_complement
from operator import itemgetter
from docopt import docopt
from hashlib import md5
from pprint import pprint

logging.basicConfig()
log = logging.getLogger(__name__)

# GTF/GFF is a standard format for DNA sequence features:
# http://www.ensembl.org/info/website/upload/gff.html
#
# The format is TSV with 9 fields as follows:
#
gtf_fields = [
    'seqname',            # e.g. chr1
    'source',             # e.g. ensembl_havana
    'feature',            # e.g. transcript, exon, CDS
    'start',              # one-based, inclusive 
    'end',                # one-based, inclusive
    'score',              # a score (unused in Ensembl GTF)
    'strand',             # + (forward), - (reverse)
    'frame',              # 0, 1, 2 
    'attribute'           # key1, val1; key2, val2; ... etc.
]
#
# The JSON template file specifies the features, fields and 
# feature-specific attributes to be extracted from the GTF 
# file, and the aliases to be applied. Example features are: 
# gene, CDS, stop_codon. Example feature-specific attributes 
# are (respectively) gene_biotype, protein_id, exon_number. 
# Fields include the annotation source and the coordinates. 
# The genome configuration CSV specifes aliases for the 
# chromosome names used in the GTF file.
#
# The are two stages to the parsing. 
#
# In the first stage a raw parse is made of the GTF file,
# producing a <feature>.csv file for each feature. These files
# apply the aliases specified in the JSON template and
# transform the coordinates to in-house zero-based [,) standard.
# The 'unique' field of the JSON template specifies a list 
# of parsed values that are expected to occur no more than 
# once in the file. An exception occurs if this is violated.
#
# The second stage (optional, but the default) will be run if 
# the JSON template specifies a minimum of gene, transcript 
# and exon and CDS features. In this stage, the output of the 
# first stage is transformed into a set of CSV files which 
# correspond closely to the database tables into which they 
# can subsequently be loaded using the appropriate manifest 
# loader. These files are therefore an external representation 
# of the core genome annotation tables in the database. 
#
# Logically the second stage could be in a separate script,
# but practically the two stages will almost always be run
# in sequence, so the two are combined here.

def main(args=None):
    """ Obtain from docopt if no args passed """

    # Parse args
    if not args:
        args = docopt(__doc__, version='GTF Parser')

    # Set logging level
    loglevel = getattr(logging, args['--log'].upper())
    log.setLevel(loglevel)

    # Log args
    log.info(str(args))

    # File paths
    input_file = os.path.abspath(args['--input'])
    template_file = os.path.abspath(args['--template'])
    config_file = os.path.abspath(args['--config'])
    work_dir = os.path.abspath(args['--output'])

    # Flags
    transform = int(args['--transform'])
    require_ccds = int(args['--require-ccds'])
    translate = int(args['--translate'])
    namespace = args['--namespace']

    # Read template
    template_dict = None
    with open(template_file) as template_fh:
        template_str = template_fh.read()
        template_dict = json.JSONDecoder(object_pairs_hook=collections.OrderedDict).decode(template_str)

    # Extract chromosome name map from the genome config file
    seqname_dict = collections.OrderedDict()
    with open(config_file) as config_fh:
        reader = csv.DictReader(config_fh)
        for config_dict in reader:
            seqname_dict[config_dict['seqname']] = config_dict['chromosome_name']

    # Stage one parse
    raw_parse(work_dir, template_dict, seqname_dict, input_file)

    # Proceed to second stage?
    if not transform:
        log.info('Second stage transform disabled, so all done')
        return        
    for feature in 'gene', 'transcript', 'exon', 'CDS':
        if not feature in template_dict:
            log.info('Feature %s required for second stage has not been parsed, quitting')
            return

    # File names for output
    db_dict = { key: "{0}/db.{1}.csv".format(work_dir, key) for key in [
        'gene', 'cds', 'cdsregion', 'transcript', 'exon']}
    
    # Stage two parse and output
    db_transform(template_dict, db_dict, require_ccds, translate, namespace)
                  

def raw_parse(work_dir, template_dict, seqname_dict, input_file):
    """ Parse specified features, fields, attributes to CSV"""

    log.info('===== Raw parse =====')

    # Features of interest
    interest = set()
    for key, obj in template_dict.items():
        interest.add(key)

    # Open CSV files, write header
    for key, template in template_dict.items():
        template['output_file'] = "{0}/{1}.csv".format(work_dir, key)
        template['output_fh'] = open(template['output_file'], 'w')
        template['writer'] = csv.DictWriter(
            template['output_fh'],
            template['attributes'].values() + template['fields'].values())
        template['writer'].writeheader()
        template['written'] = 0

    # If a unique key tuple has been defined, prepare a dejavu set
    for key, template in template_dict.items():
        if 'unique' in template:
            template['dejavu'] = set()
        
    # Parse for features of interest
    log.info('Reading %s...', input_file)
    with gzip.open(input_file) as input_fh:
        # Uncertain number of header lines, so awkward to use csv.reader
        for line in input_fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) != len(gtf_fields):
                raise Exception("Unexpected number of columns (%d not %d) in gtf file line %s\n", 
                                % (len(parts), len(gtf_fields), line))
            feature_dict = dict(zip(gtf_fields, parts))
            # Handle specified features for configured chromosomes/sequences, ignore otherwise
            if (feature_dict['feature'] in interest and
                feature_dict['seqname'] in seqname_dict):
                _handle(template_dict, feature_dict, seqname_dict, line)

    # Output
    for key, template in template_dict.items():
        template['output_fh'].close()
        log.info('Wrote %d items to %s', 
                 template['written'], template['output_file'])


def _handle(template_dict, feature_dict, seqname_dict, line_for_log):
    """
    Instantiate a new dict collecting the templated fields and 
    attributes from the feature_dict, applying aliases and 
    transforming the coordinates. Write to file.
    """
    
    # Unpack the attribute field
    attribute_regex = r"(.*?)\s+\"(.*?)\";\s?"
    attribute_dict = {}
    if 'attribute' in feature_dict:
        kv_pairs = re.findall(attribute_regex, feature_dict['attribute'])
        attribute_dict = collections.OrderedDict(kv_pairs)

    # Do we have a template for this feature?
    feature = feature_dict['feature']
    if not feature in template_dict:
        raise Exception('No template provided for %s' % (feature))
    template = template_dict[feature]

    # The target
    templated_dict = collections.OrderedDict()

    # Transform coordinates in place, overwriting feature_dict (..)
    _transform_coordinates(feature_dict, seqname_dict, line_for_log)
    
    # Collect specified fields and attributes, apply aliases
    for source_dict, alias_dict in [
            (attribute_dict, template['attributes']),
            (feature_dict, template['fields'])
    ]:
        for from_key in sorted(alias_dict.keys()):
            if from_key in source_dict:
                to_key = alias_dict[from_key]
                if source_dict[from_key] != '.':
                    templated_dict[to_key] = source_dict[from_key]

    # If unique keys not specified, no check to do - just write and return
    if not 'unique' in template:
        template['writer'].writerow(templated_dict)
        template['written'] += 1
        return
        
    # Check the tuple of fields that are meant to be unique
    tup = tuple(map(templated_dict.get, template['unique']))
    if tup in template['dejavu']:
        raise Exception("Error %s %s occurs more than once\n%s" % 
                        (str(template['unique']), str(tup), line_for_log))

    # Passes, record tuple and write out
    template['dejavu'].add(tup)
    template['writer'].writerow(templated_dict)
    template['written'] += 1


def _hash_dict(input_dict):
    inst = md5()
    inst.update(str(input_dict))
    return inst.hexdigest()


def db_transform(template_dict, db_dict, require_ccds, translate, namespace):
    """
    Transform data from raw parse to produce CSV files for database load
    """

    log.info('===== Reading data for transform =====')

    # Read gene file
    gene_dict = {}
    with open(template_dict['gene']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for row_dict in reader:
            gene_dict[row_dict['gene_accession']] = row_dict
    log.info('Loaded %d genes', len(gene_dict))

    # Add a namespace for gene if not provided
    for gene_accession, gene in gene_dict.items():
        if not 'namespace' in gene or not gene['namespace']:
            gene['namespace'] = namespace

    # Read transcripts file
    transcript_dict = {}
    with open(template_dict['transcript']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for row_dict in reader:
            transcript_dict[row_dict['transcript_accession']] = row_dict
    log.info('Loaded %d transcripts', len(transcript_dict))

    # Read the coding exons into a handy dict of dicts
    # To access exon_dict for transcript_order within transcript 
    # transcript_accession use
    # exon_dict = cds_dict[transcript_accesson][transcript_order]
    cds_dict = collections.OrderedDict()
    with open(template_dict['CDS']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for exon_dict in reader:
            transcript_accession = exon_dict['transcript_accession']
            transcript_order = int(exon_dict['transcript_order'])
            if not transcript_accession in cds_dict:
                cds_dict[transcript_accession] = collections.OrderedDict()
            cds_dict[transcript_accession][transcript_order] = exon_dict
        input_fh.close()
    log.info('Loaded %d CDS', len(cds_dict))

    # Read stop codons if present
    # This is a dict of codon_dicts keyed on tuple 
    # (transcript_accession, transcript_order)
    stop_dict = {}
    if 'stop_codon' in template_dict:
        with open(template_dict['stop_codon']['output_file'], 'r') as input_fh:
            reader = csv.DictReader(input_fh)
            for codon_dict in reader:
                transcript_accession = codon_dict['transcript_accession']
                transcript_order = int(codon_dict['transcript_order'])
                tup = transcript_accession, transcript_order
                if tup in stop_dict:
                    raise Exception('Multiple stop features for same codon %s', str(tup))
                stop_dict[tup] = codon_dict
            input_fh.close()
        log.info('Loaded %d stop codons', len(stop_dict))

    log.info('===== Cross feature checks =====')
                                 
    # Is there a gene_accession for every transcript?
    zap = set(gene_dict.keys())
    for transcript_accession, transcript in transcript_dict.items():
        if not 'gene_accession' in transcript:
            raise Exception('No gene accession provided for %s', transcript_accession)
        gene_accession = transcript['gene_accession'] 
        if not gene_accession in gene_dict:
            raise Exception('No gene feature was loaded for gene %s referred to in transcript %s', 
                            gene_accession, transcript_accession)
        zap.discard(gene_accession)
        gene = gene_dict[gene_accession]
        # Are the transcript coordinates within the gene coordinates?
        check_contained('Transcript', transcript, transcript_accession,
                        'gene', gene, gene_accession)
    # Did we find a transcript for all genes?
    if len(zap):
        raise Exception('Gene %s has no associated transcript' % (zap.pop()))
        
    # Did we read a transcript feature for every coding exon?
    for transcript_accession, cds in cds_dict.items():
        for transcript_order, exon_dict in cds.items():
            if not transcript_accession in transcript_dict:
                 raise Exception('No transcript feature read for CDS exon %s:%s' %
                                 (transcript_accession, transcript_order))
            # Are the exon coordinates within the transcript coordinates?
            transcript = transcript_dict[transcript_accession]
            if not 'transcript_order' in exon_dict:
                raise Exception('No exon number in %s\n' % (str(exon_dict)))
            check_contained('CDS exon', exon_dict, exon_dict['transcript_order'],
                            'transcript', transcript, transcript_accession)

    # Did we find a transcript and CDS for every stop codon?
    for (transcript_accession, transcript_order), codon_dict in stop_dict.items():
        if not transcript_accession in cds_dict:
            raise Exception('Failed to find transcript for stop codon %s:%d' % 
                            (transcript_accession, codon_dict['transcript_order']))
            # Are the stop codon coordinates within the transcript coordinates?
            transcript = transcript_dict[transcript_accession]
            check_contained('Stop codon', codon_dict, codon_dict['transcript_order'],
                            'transcript', transcript, transcript_accession)

    # Did we read a transcript feature for every exon?
    with open(template_dict['exon']['output_file'], 'r') as input_fh:
        # Exon data not read in since no real need and can be large
        reader = csv.DictReader(input_fh)
        for exon_dict in reader:
            transcript_accession = exon_dict['transcript_accession']
            if not transcript_accession in transcript_dict:
                 raise Exception('No transcript feature read for exon %s:%s' %
                                 (transcript_accession, transcript_order))
            # Are the exon coordinates within the transcript coordinates?
            transcript = transcript_dict[transcript_accession]
            check_contained('Exon', exon_dict, exon_dict['transcript_order'],
                            'transcript', transcript, transcript_accession)

    # Overlap/ordering checks
    for transcript_accession, cds in cds_dict.items():
        # Order by exon_number
        order = sorted(cds.keys())
        if int(transcript_dict[transcript_accession]['strand']) == -1:
            order = sorted(order, reverse=True)
        # Are the coding exons in transcript_order in the sequence?
        prev_dict = None
        for transcript_order in order:
            exon_dict = cds[transcript_order]
            if not prev_dict:
                prev_dict = exon_dict
                continue
            if int(exon_dict['start']) <= int(prev_dict['start']):
                print transcript
                raise Exception('Transcript/sequence order inconsistency for CDS exons %s %s at\n   %s\nand %s %s at\n   %s\n' %
                                (transcript_accession, exon_dict['transcript_order'], str(exon_dict),
                                 transcript_accession, prev_dict['transcript_order'], str(prev_dict)))
        # Are there any overlaps?
        prev_dict = None
        for transcript_order in order:
            exon_dict = cds[transcript_order]
            if not prev_dict:
                prev_dict = exon_dict
                continue
            if int(exon_dict['start']) < int(prev_dict['end']):
                print transcript
                raise Exception('Overlapping CDS exons %s %s at\n   %s\nand %s %s at\n   %s\n' %
                                (transcript_accession, exon_dict['transcript_order'], str(exon_dict),
                                 transcript_accession, prev_dict['transcript_order'], str(prev_dict)))

    log.info('===== Applying transformations =====')

    # Integrate the stop codons
    # Important to do this before calculating phases
    for (transcript_accession, transcript_order), codon_dict in stop_dict.items():
        if not transcript_accession in cds_dict:
            raise Exception('Failed to find %s', transcript_accession)
        if transcript_order in cds_dict[transcript_accession]:
            # Append to the existing codon
            # Recast to str for consistency with other CDS and stop_codons
            exon_dict = cds_dict[transcript_accession][transcript_order]
            exon_dict['start'] = str(min(int(exon_dict['start']),
                                         int(codon_dict['start'])))
            exon_dict['end'] = str(max(int(exon_dict['end']),
                                       int(codon_dict['end'])))
        else:
            # Append to list of codons for this transcript
            log.warn('Appending exon with only stop codon to %s', transcript_accession)
            cds_dict[transcript_accession][transcript_order] = codon_dict

    # Assign cds_order
    # We could use the transcript_order for ordering, but that
    # doesn't always start from 1
    for transcript_accession, cds in cds_dict.items():
        order = 1
        for transcript_order in sorted(cds.keys()):
            cds[transcript_order]['cds_order'] = order
            order += 1

    # Calculate phases
    # Frame is the offset of the 1st codon from the feature start:
    # http://www.ensembl.org/info/website/upload/gff.html
    # Phase is the offset of the feature start from the 1st codon:
    # https://github.com/Ensembl/ensembl/blob/release/79/modules/Bio/EnsEMBL/Exon.pm
    # Hence phase = ( - frame ) % 3
    # See also glossary
    # http://www.ensembl.org/Help/Glossary
    for transcript_accession, cds in cds_dict.items():
        for transcript_order, exon_dict in cds.items():
            exon_length = int(exon_dict['end']) - int(exon_dict['start'])
            exon_dict['phase_start'] = (- int(exon_dict['frame'])) % 3
            exon_dict['phase_end'] = (exon_dict['phase_start'] + exon_length) % 3

    if translate:
        raise Exception('Translate option not yet implemented')
        pass
                       
    # Calculate cds and cdsregions data
    #
    # The cds and cdsregion tables in the database represent
    # coding sequences i.e. the parts of a genomic sequence
    # that code for a protein. This includes the start and stop
    # codons. The cdsregion table holds the coordinates of the
    # coding regions within individual exons. The cds table 
    # is the transcript-level grouping of cdsregions.
    #
    # Because transcripts have non-coding exons and/or exons 
    # which are partly coding and partly non-coding, multiple 
    # transcripts can have the same coding sequence.
    #
    # In the following section we generate cds and cdsregion
    # records via a process of coordinate de-duplication.
    # We say that transcripts have the same coding sequence 
    # if their coding regions are on the same strand, in the 
    # same order and with the same coordinates. The script
    # identifies the distinct coding sequences and generates
    # the one-to-many mappings from coding sequence to 
    # transcript.
    #
    # Separately, the CCDS project 
    # http://www.ensembl.org/info/genome/genebuild/ccds.html
    # provides a curated set of coding sequences, for human
    # and mouse. Each coding sequence in the CCDS set carries
    # a CCDS id and the coordinates of the coding regions.
    #
    # The Ensembl GTF files for human and mouse include the 
    # CCDS id as an attribute of the transcript feature. 
    # The naive expectation would be that in those GTF files,
    # transcripts with the same CCDS id should have the same
    # coding sequence, according to our definition above. 
    # In fact, this is almost always the case. The handful of
    # exceptions are logged. One might also expect that 
    # every coding transcript in the human and mouse GTF files 
    # should have a CCDS id, but this is not the case at all. 
    # There are other 'quality criteria' used for the 
    # assignment of CCDS ids and about 1000 genes in the Ensembl
    # set have no transcript with a CCDS id.
    # 
    # Transcripts that have a coding sequence but no CCDS id
    # (e.g. not human or mouse, or human/mouse but no CCDS id 
    # assigned) are handled in this script according to the 
    # '--require-ccds' flag. The flag should be set false for
    # genomes such as the rat genome for which there is no CCDS
    # data. This will ensure that data fields required for the 
    # front end are populated with appropriate alternative 
    # information.

    # Maps to enable CCDS handling
    ccds_to_transcript = {}
    transcript_to_ccds = {}
    for transcript_accession, transcript in transcript_dict.items():
        if 'ccds_accession' in transcript and transcript['ccds_accession']:
            ccds_accession = transcript['ccds_accession']
            if not ccds_accession in ccds_to_transcript:
                ccds_to_transcript[ccds_accession] = []
            ccds_to_transcript[ccds_accession].append(transcript_accession)
            transcript_to_ccds[transcript_accession] = ccds_accession

    # Maps between coding sequences and transcripts
    cds_to_transcript = {}
    transcript_to_cds = {}
    count = 0
    for transcript_accession, cds in cds_dict.items():
        if require_ccds:
            # Ignore transcripts without a ccds_id, don't even
            # calculate the coding sequence
            if not transcript_accession in transcript_to_ccds:
                continue
        # The following hash defines uniqueness for the cds file. It
        # depends on both coordinates and other fields, so there is
        # possible non-uniqueness at the level of coordinates alone.
        # In practise, almost all records will be unique in terms of
        # coordinates.
        dict_to_hash = collections.OrderedDict({
            'name': None,
            'namespace': None,
            'is_consensus': 0,
            'exon_coords': []
        })
        # Append coords
        for transcript_order in sorted(cds.keys()):
            exon_dict = cds[transcript_order]
            coord_dict = collections.OrderedDict({key: exon_dict[key] for key in (
                'chromosome_name', 'start', 'end', 'strand',
                'phase_start', 'phase_end',
                'cds_order')})
            dict_to_hash['exon_coords'].append(coord_dict)
        # Attached CCDS id if available
        if not dict_to_hash['name'] and transcript_accession in transcript_to_ccds:
            dict_to_hash['name'] = transcript_to_ccds[transcript_accession]
            dict_to_hash['namespace'] = 'ccds'
            dict_to_hash['is_consensus'] = 1
        # Else attach a protein id
        if not dict_to_hash['name'] and transcript_accession in cds_dict:
            cds = cds_dict[transcript_accession]
            for transcript_order, exon_dict in cds.items():
                if 'protein_id' in exon_dict:
                    dict_to_hash['name'] = exon_dict['protein_id']
                    break
        # Else create an id
        if not dict_to_hash['name']:
            dict_to_hash['name'] = 'CDS%06d' % (len(cds_to_transcript) + 1)
        hash_val = _hash_dict(dict_to_hash)
        # Have we got it?
        if not hash_val in cds_to_transcript:
            # No, add it
            cds_to_transcript[hash_val] = {
                'cds_accession': "eq_{0}".format(len(cds_to_transcript)), 
                'namespace': dict_to_hash['namespace'],
                'name': dict_to_hash['name'],
                'gene_accession': transcript_dict[transcript_accession]['gene_accession'],
                'is_consensus': dict_to_hash['is_consensus'],
                'status': None,
                'transcript_accession_list': [],
                'exon_coords': dict_to_hash['exon_coords']
            }
            count += len(dict_to_hash['exon_coords'])
        # Cross reference cds/transcript
        cds_to_transcript[hash_val]['transcript_accession_list'].append(transcript_accession)
        transcript_to_cds[transcript_accession] = hash_val
    log.info('Found %d coding sequences with total %d exons', 
              len(cds_to_transcript), count)

    # Calculate cds coordinates as the envelope of it's constituent exons
    for hash_val, obj in cds_to_transcript.items():
        obj['chromosome_name'] = obj['exon_coords'][0]['chromosome_name']
        obj['strand'] = obj['exon_coords'][0]['strand']
        start = None
        end = None
        for coord_dict in obj['exon_coords']:
            if not start or int(coord_dict['start']) < start:
                start = int(coord_dict['start'])
            if not end or int(coord_dict['end']) > end:
                end = int(coord_dict['end'])
        # Cast to string for consistency
        obj['start'] = str(start)
        obj['end'] = str(end)

    # For logging purposes only
    # Check that transcripts with same CCDS have the same coding sequence
    for ccds_accession, transcript_accession_list in ccds_to_transcript.items():
        first_transcript_accession = transcript_accession_list[0]
        first_hash_val = transcript_to_cds[first_transcript_accession]
        for transcript_accession in transcript_accession_list:
            if not transcript_to_cds[transcript_accession] == first_hash_val:
                log.warning('Transcripts %s and %s with same CCDS have different coding sequences',
                            first_transcript_accession,
                            transcript_accession)

    # Assign namespace to any transcripts without a namespace
    for transcript_accession, transcript in transcript_dict.items():
        if 'namespace' not in transcript:
            transcript['namespace'] = namespace

    # Write the cds name into the ccds_accession field of transcript,
    # if that field is not set. This is a bit of an abuse but the 
    # field is needed for front end display.
    for transcript_accession, transcript in transcript_dict.items():
        if not 'ccds_accession' in transcript or not transcript['ccds_accession']:
            if transcript_accession in transcript_to_cds:
                transcript['ccds_accession'] = cds_to_transcript[transcript_to_cds[transcript_accession]]['name']

    # Label each transcript with it's cds accession if available.
    for transcript_accession, transcript in transcript_dict.items():
        transcript['cds_accession'] = None
        if transcript_accession in transcript_to_cds:
            transcript['cds_accession'] = cds_to_transcript[transcript_to_cds[transcript_accession]]['cds_accession']

    log.info('===== Output =====')

    # Produce the gene file
    fieldnames = [
        'gene_accession',
        'name',
        'namespace',
        'biotype',
        'chromosome_name',
        'start',
        'end',
        'strand'
    ]
    with open(db_dict['gene'], 'w') as output_fh:
        writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()
        count = 0
        for gene_accession, gene in gene_dict.items():
            row_dict = {key: gene[key] for key in fieldnames}
            writer.writerow(row_dict)
            count += 1
        output_fh.close()
        log.info('Wrote %d regions to %s', count, db_dict['gene'])

    # Produce the cds file
    fieldnames = [
        'cds_accession',
        'name',
        'namespace',
        'gene_accession',
        'is_consensus',
        'status',
        'chromosome_name',
        'start',
        'end',
        'strand'
    ]
    with open(db_dict['cds'], 'w') as output_fh:
        writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()
        count = 0
        for hash_val, obj in cds_to_transcript.items():
            row_dict = {key: obj[key] for key in fieldnames}
            writer.writerow(row_dict)
            count += 1
        output_fh.close()
        log.info('Wrote %d regions to %s', count, db_dict['cds'])

    # Produce the cdsregions file
    fieldnames = [
        'cds_accession',
        'cds_order',
        'phase_start',
        'phase_end',
        'chromosome_name',
        'start',
        'end',
        'strand'
   ]
    with open(db_dict['cdsregion'], 'w') as output_fh:
        writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()
        count = 0
        for hash_val, obj in cds_to_transcript.items():
            for coord_dict in obj['exon_coords']:
                coord_dict['cds_accession'] = obj['cds_accession']
                row_dict = {key: coord_dict[key] for key in fieldnames}
                writer.writerow(row_dict)
                count += 1
        output_fh.close()
        log.info('Wrote %d exons to %s', count, db_dict['cdsregion'])

    # Produce the transcripts file
    fieldnames = [
        'transcript_accession',
        'gene_accession',
        'ccds_accession',
        'name',
        'namespace',
        'biotype',
        'chromosome_name',
        'start',
        'end',
        'strand',
        'transcript_support_level',
        'cds_accession'
    ]
    with open(db_dict['transcript'], 'w') as output_fh:
        writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()
        count = 0
        for transcript_accession, transcript in transcript_dict.items():
            row_dict = {key: transcript[key] for key in fieldnames}
            writer.writerow(row_dict)
            count += 1
        output_fh.close()
        log.info('Wrote %d regions to %s', count, db_dict['transcript'])

    # Exon data are not read in since no real need and can be large.
    # No transformations, but extra columns are included when the
    # exon has a coding region.
    # REGRESSION NOTE
    # In the current database, if the transcript has a corresponding
    # entry in the cds table, the exon phase start/ends are OVERWRITTEN 
    # with the coding coordinate start/end. The column 
    # 'phase_db_overwrite' in the output file indicates to the loader 
    # when to do this.
    fieldnames = [
        'transcript_accession',
        'transcript_order',
        'phase_start',
        'phase_end',
        'chromosome_name',
        'start',
        'end',
        'strand',
        'coding_start',
        'coding_end',
        'phase_db_overwrite'
    ]
    with open(db_dict['exon'], 'w') as output_fh:
        writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()
        count = 0
        with open(template_dict['exon']['output_file'], 'r') as input_fh:
            reader = csv.DictReader(input_fh)
            for exon_dict in reader:
                transcript_accession = exon_dict['transcript_accession']
                transcript_order = int(exon_dict['transcript_order'])
                exon_dict['phase_start'] = -1
                exon_dict['phase_end'] = -1
                exon_dict['coding_start'] = -1
                exon_dict['coding_end'] = -1
                exon_dict['phase_db_overwrite'] = 0
                if transcript_accession in cds_dict:
                    cds = cds_dict[transcript_accession]
                    if transcript_order in cds:
                        # This exon contains a coding region
                        coding_exon_dict = cds[transcript_order]
                        # If the exon is coding at an end assign the coding phase
                        if int(exon_dict['strand']) == 1:
                            if int(exon_dict['start']) == int(coding_exon_dict['start']):
                                exon_dict['phase_start'] = coding_exon_dict['phase_start']
                            if int(exon_dict['end']) == int(coding_exon_dict['end']):
                                exon_dict['phase_end'] = coding_exon_dict['phase_end']
                        else:
                            if int(exon_dict['start']) == int(coding_exon_dict['start']):
                                exon_dict['phase_end'] = coding_exon_dict['phase_end']
                            if int(exon_dict['end']) == int(coding_exon_dict['end']):
                                exon_dict['phase_start'] = coding_exon_dict['phase_start']
                        # The end of the last exon is non-coding so the phase should be -1
                        # http://www.ensembl.org/Help/Glossary
                        if int(coding_exon_dict['cds_order']) == len(cds):
                            exon_dict['phase_end'] = -1
                        exon_dict['coding_start'] = coding_exon_dict['start']
                        exon_dict['coding_end'] = coding_exon_dict['end']
                        # Loader overwrite flag for transcripts with a cds
                        if transcript_accession in transcript_to_cds:
                            exon_dict['phase_db_overwrite'] = 1
                if int(exon_dict['start']) > int(exon_dict['end']):
                    log.warn('Omitting invalid exon coordinates for transcript %s', transcript_accession)
                    continue
                writer.writerow(exon_dict)
                count += 1
        output_fh.close()
        log.info('Wrote %d exons to %s', count, db_dict['exon'])


def check_contained(feature_one, coords_one, accession_one,
                    feature_two, coords_two, accession_two):
    # Raise exception if coords_one are not within coords two
    if (coords_one['chromosome_name'] != coords_two['chromosome_name'] or
        coords_one['strand'] != coords_two['strand'] or
        int(coords_one['start']) < int(coords_two['start']) or
        int(coords_one['end']) > int(coords_two['end'])):
        raise Exception('%s %s at\n   %s\nis not contained within the coordinates of %s %s at\n   %s' % 
                        (feature_one, accession_one, str(coords_one),
                         feature_two, accession_two, str(coords_two)))

def _transform_coordinates(feature_dict, seqname_dict, line_for_log):

    # Modifications in place, input overwritten

    # Sequence name must be in the config
    if not feature_dict['seqname'] in seqname_dict:
        raise Exception('Invalid seqname, not in config\n%s' % (line_for_log))
    feature_dict['seqname'] = seqname_dict[feature_dict['seqname']]

    # Strand must be '+' or '-'
    if not feature_dict['strand'] in ('+', '-'):
        raise Exception('Invalid strand\n%s' % (line_for_log))
    feature_dict['strand'] = 1 if feature_dict['strand'] == '+' else -1

    # Validate coordinates
    # Checks uses pre-transformed, one-based [,] per GTF specification
    if int(feature_dict['end']) < int(feature_dict['start']):
        raise Exception('Inverted coordinates\n%s' % (line_for_log))
    if int(feature_dict['start']) < 1:
        raise Exception('Start coordinate out of bounds\n%s' % (line_for_log))
    # Transform to zero-based [,) and recast to str for consistency
    feature_dict['start'] = str(int(feature_dict['start']) - 1)


if __name__ == '__main__':
    main()

