"""Parse GTF for Gene, Transcript, Exon and CDS features

Usage:
  parse_gtf.py (--input=<file.gz>) (--template=<file>) (--config=<file>) [--check=<bool>] [--withstop=<bool>] [--log=<loglevel>]

Options:
  --help                    Show this screen
  --input=<name>            GTF file, gzipped
  --template=<file>         JSON template for I/O field mapping
  --config=<file>           CSV file restricting seqnames to extract features for
  --check=<bool>            Check redundant data is identical [default: 0]
  --withstop=<bool>         Merge stop_codons with CDS exons [default: 1]
  --log=<loglevel>          Logging level [default: ERROR]

"""
import os
import re
import csv
import gzip
import logging
import collections
import json

from operator import itemgetter

from docopt import docopt
from hashlib import md5
from pprint import pprint

logging.basicConfig()
log = logging.getLogger(__name__)

def _hash_dict(input_dict):
    inst = md5()
    inst.update(str(input_dict))
    return inst.hexdigest()
    

def _handle(template, check, feature_dict):

    # Unpack the attribute field 
    attribute_regex = r"(.*?)\s+\"(.*?)\";\s?"
    attribute_dict = {}
    if 'attribute' in feature_dict:
        kv_pairs = re.findall(attribute_regex, feature_dict['attribute'])
        attribute_dict = collections.OrderedDict(kv_pairs)

    feature = feature_dict['feature']
    
    for gtf_key, obj in template.items():
        if feature == gtf_key:

            model = collections.OrderedDict()

            for source_dict, map_dict in [
                    (attribute_dict, obj['attributes']),
                    (feature_dict, obj['fields'])
            ]:
                for from_key in sorted(map_dict.keys()):
                    if from_key in source_dict:
                        to_key = map_dict[from_key]
                        if source_dict[from_key] != '.':
                            model[to_key] = source_dict[from_key]

            if not 'unique' in obj:
                obj['writer'].writerow(model)
                continue

            check_value = None
            if check:
                check_value = _hash_dict(model)
            tup = tuple(map(model.get, obj['unique']))
            if tup in obj['dejavu']:
                if 'enforce_unique' in obj and obj['enforce_unique']:
                    raise Exception("Error %s occurs more than once" % str(tup))
                else:
                    if check:
                        assert check_value == obj['dejavu'][tup]
            else:
                obj['dejavu'][tup] = check_value
                obj['writer'].writerow(model)


def parse(template_file, config_file, check, input_file):
    
    # GTF format is a TSV format file containing gene features. 
    # Reference:
    # http://www.ensembl.org/info/website/upload/gff.html

    # The fields are    
    gtf_fields = [
        'seqname',            # e.g. chromosome_name
        'source',             # e.g. ensembl_havana
        'feature',            # feature name e.g. CDS, transcript, exon
        'start',              
        'end',
        'score',
        'strand',             # + (forward), - (reverse)
        'frame',
        'attribute'           # key1, val1; key2, val2; ... etc.
    ]

    # Read template
    template = None
    with open(template_file) as template_fh:
        template_str = template_fh.read()
        template =json.JSONDecoder(object_pairs_hook=collections.OrderedDict).decode(template_str)

    # Read template
    seqname_dict = {}
    with open(config_file) as config_fh:
        reader = csv.DictReader(config_fh)
        for config_dict in reader:
            seqname_dict[config_dict['seqname']] = config_dict['chromosome_name']
            
    # Features of interest
    interest = set()
    for key, obj in template.items():
        interest.add(key)

    # Open CSV files, write header
    for key, obj in template.items():
        obj['output_file'] = "{0}.csv".format(key)
        obj['output_fh'] = open(obj['output_file'], 'w')
        obj['writer'] = csv.DictWriter(
            obj['output_fh'], 
            obj['attributes'].values() + obj['fields'].values())
        obj['writer'].writeheader()

    # If a unique field name tuple has been defined, prepare a dejavu dict
    for key, obj in template.items():
        if 'unique' in obj:
            obj['dejavu'] = {}

    # Parse for features of interest
    log.info('Reading %s...', input_file)
    with gzip.open(input_file) as input_fh:
        # Uncertain number of header lines, tricky to use csv.reader
        for line in input_fh:
            if line.startswith('#'):
                continue
            stripped = line.rstrip()
            feature_dict = dict(zip(gtf_fields, stripped.split('\t')))
            if (feature_dict['feature'] in interest and
                feature_dict['seqname'] in seqname_dict):
                _transform_coordinates(feature_dict, seqname_dict)
                _handle(template, check, feature_dict)
#                if len(template['gene']['dejavu']) > 100:
#                    break
                
    # Output
    for key, obj in template.items():
        obj['output_fh'].close()
        log.info('Wrote %d items to %s', len(obj['dejavu']), obj['output_file'])


def _transform_coordinates(feature_dict, seqname_dict):
    feature_dict['seqname'] = seqname_dict[feature_dict['seqname']]
    feature_dict['strand'] = 1 if feature_dict['strand'] == '+' else -1
    feature_dict['start'] = str(int(feature_dict['start']) - 1)


def dbmap(transcript_file, exon_file, CDS_file, stop_file, 
          db_exon_file, db_cdsregion_file):
    """
    Take the CSV files of parsed GTF features, transform into 
    CSV files corresponding directly to the database tables.
    """

    # Transcripts with CCDS
    with_ccds = set()
    with open(transcript_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for transcript_dict in reader:
            if 'ccds_accession' in transcript_dict:
                with_ccds.add(transcript_dict['accession'])
    
    # Read stop codons
    stop_dict = {}
    with open(stop_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for codon_dict in reader:
            tup = codon_dict['transcript_name'], codon_dict['transcript_order']
            if tup in stop_dict:
                raise Exception('Multiple stop codons for %s' % 
                                (", ".join(tup)))
            stop_dict[tup] = codon_dict
        input_fh.close()
    log.debug('Loaded %d stop codons', len(stop_dict))

    # Read coding exons
    # Grab fields names
    cds_dict = collections.OrderedDict()
    fieldnames = None
    with open(CDS_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for exon_dict in reader:
            if not fieldnames:
                fieldnames = reader.fieldnames
            transcript_name = exon_dict['transcript_name']
            transcript_order = exon_dict['transcript_order']
            if not transcript_name in cds_dict:
                cds_dict[transcript_name] = collections.OrderedDict()
            cds_dict[transcript_name][int(transcript_order)] = exon_dict
        input_fh.close()
    log.debug('Loaded %d CDS exons', len(cds_dict))

    # Assign cds_order
    # A bit unecessary, could use transcript_order directly
    for transcript_name, cds in cds_dict.items():
        order = 1
        for transcript_order in sorted(cds.keys()):
            cds[transcript_order]['cds_order'] = order
            order += 1

    # Append stop codons and calculate phases 
    for transcript_name, cds in cds_dict.items():
        for transcript_order, exon_dict in cds.items():
            tup = transcript_name, transcript_order
            # Stop codon
            if tup in stop_dict:
                codon_dict = stop_dict[tup]
                exon_dict['start'] = min(int(exon_dict['start']), 
                                         int(codon_dict['start']))
                exon_dict['end'] = max(int(exon_dict['end']), 
                                       int(codon_dict['end']))
            # Insert phases
            exon_length = int(exon_dict['end']) - int(exon_dict['start'])
            exon_dict['phase_start'] = ( 3 - int(exon_dict['frame']) ) % 3
            exon_dict['phase_end'] = ( exon_dict['phase_start'] + exon_length ) % 3

    # Fix the fields
    fieldnames.append('cds_order')
    fieldnames.append('phase_start')
    fieldnames.append('phase_end')

    # Write cdsregions
    with open(db_cdsregion_file, 'w') as output_fh:
        writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
        writer.writeheader()
        for transcript_name, cds in cds_dict.items():
            # Only write exons for the CCDS transcripts
            if transcript_name in with_ccds:
                for transcript_order, exon_dict in cds.items():
                    writer.writerow(exon_dict)
        output_fh.close() 

    fieldnames = None
    with open(db_exon_file, 'w') as output_fh:
        writer = None
        with open(exon_file, 'r') as input_fh:
            reader = csv.DictReader(input_fh)
            for exon_dict in reader:
                # Insert phases
                exon_dict['phase_start'] = -1
                exon_dict['phase_end'] = -1
                # These fields will actually hold the coding exon 
                # start/ends if coding, not the phases!
                if exon_dict['accession'] in cds_dict:
                    cds = cds_dict[exon_dict['accession']]
                    transcript_order = int(exon_dict['transcript_order'])
                    if transcript_order in cds:
                        coding_exon_dict = cds[transcript_order]
                        exon_dict['phase_start'] = coding_exon_dict['start']  # Nasty
                        exon_dict['phase_end'] = coding_exon_dict['end']      # Nasty
                if not writer:
                    fieldnames = reader.fieldnames
                    fieldnames.append('phase_start')
                    fieldnames.append('phase_end')
                    writer = csv.DictWriter(output_fh, fieldnames=fieldnames)
                    writer.writeheader()
                writer.writerow(exon_dict)
        output_fh.close() 
        

def main():

    # Parse args
    args = docopt(__doc__, version='Breakpoints Annotator')

    # Set logging level
    loglevel = getattr(logging, args['--log'].upper())
    log.setLevel(loglevel)

    # Log args
    log.info(str(args))

    # File paths
    input_file = os.path.abspath(args['--input'])
    template_file = os.path.abspath(args['--template'])
    config_file = os.path.abspath(args['--config'])

    # Flags
    check = args['--check']
    withstop = args['--withstop']

    # Parse
    parse(template_file, config_file, check, input_file)

    # Merge stops
    if withstop:
        dbmap('transcript.csv', 'exon.csv', 'CDS.csv', 'stop_codon.csv',
              'db.exon.csv', 'db.cdsregion.csv')

if __name__ == '__main__':
    main()

