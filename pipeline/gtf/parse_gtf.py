"""Parse GTF for Gene, Transcript, Exon and CDS features

Usage:
  parse_gtf.py (--input=<file.gz>) (--template=<file>) (--config=<file>) (--namespace=<str>) [--transform=<bool>] [--require-ccds=<bool>] [--log=<loglevel>]

Options:
  --help                    Show this screen
  --input=<name>            GTF file, gzipped
  --template=<file>         JSON template for I/O field mapping
  --config=<file>           Genome config, CSV file specifying the seqnames to extract features for
  --namespace=<str>         Where the data comes from
  --transform=<bool>        Prepare files ready for database loading [default: 1]
  --require-ccds=<bool>     Output cds/cdsregions only if there is a CCDS identifier [default: 1]
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

# GFF/GTF is a standard format for DNA sequence features:
# http://www.ensembl.org/info/website/upload/gff.html
#
# The format is TSV with 9 fields as follows:
#
gtf_fields = [
    'seqname',            # e.g. chr1
    'source',             # e.g. ensembl_havana
    'feature',            # feature name e.g. CDS, transcript, exon
    'start',              # one-based [,]
    'end',                # 
    'score',              # Unused in Ensembl GFF
    'strand',             # + (forward), - (reverse)
    'frame',              # 0, 1, 2 (frame + phase) % 3 = 0
    'attribute'           # key1, val1; key2, val2; ... etc.
]
#
# The parser takes a JSON template file which specifies the 
# features, common fields and feature-specific attributes 
# to be extracted from the GTF file, and the aliases to be 
# applied. Example features are: gene, CDS, stop_codon. 
# Example feature-specific attributes are (respectively)
# gene_biotype, protein_id, exon_number. Common fields 
# include the annotation source and the coordinates. A
# second config file is a genome configuration CSV file which 
# specifes aliases for the chromosome names use in the 
# GTF file.
#
# The are two stages to the parsing. 
#
# In the first stage a raw parse is made of the GTF file,
# producing one <feature>.csv file per feature. These files
# apply the aliases specified in the JSON template and
# transform the coordinates to zero-based [,). The 'unique' 
# field of the JSON template specifies a list of 
# field/attribute values that are expected to occur no 
# more than once in the file. An exception occurs if this
# is violated.
#
# The second stage (optional) can be run if the JSON template 
# specifies a minimum of gene, transcript and exon and CDS 
# features. It is to transform the output of the first stage
# into a set of CSV files which correspond closely to the
# database tables into which they can be loaded. These files
# are therefore an external representation of the core 
# genome annotations in the database. 
#
# Logically the second stage could be in a separate script,
# but practically the two stages will almost always be run
# in sequence, so the two are combined.

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
    transform = int(args['--transform'])
    require_ccds = int(args['--require-ccds'])
    namespace = args['--namespace']

    # Read template
    template = None
    with open(template_file) as template_fh:
        template_str = template_fh.read()
        template = json.JSONDecoder(object_pairs_hook=collections.OrderedDict).decode(template_str)

    # Chromosome map using genome config file
    seqname_dict = collections.OrderedDict()
    with open(config_file) as config_fh:
        reader = csv.DictReader(config_fh)
        for config_dict in reader:
            seqname_dict[config_dict['seqname']] = config_dict['chromosome_name']

    # Stage one parse
    raw_parse(template, seqname_dict, input_file)

    # Proceed to second stage?
    if not transform:
        log.info('Transform disabled, returning')
        return        
    for feature in 'gene', 'transcript', 'exon', 'CDS':
        if not feature in template:
            log.info('Feature %s not parsed, returning')
            return

    # File names for output
    db_dict = { key: "db.{0}.csv".format(key) for key in [
        'gene', 'cds', 'cdsregion', 'transcript', 'exon']}
    
    # Stage two parse
    db_transform(template, db_dict, require_ccds, namespace)
                  

def raw_parse(template, seqname_dict, input_file):
    """ Parse specified features, fields, attributes to CSV"""

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

    # If a unique field name tuple has been defined, prepare a dejavu set
    for key, obj in template.items():
        if 'unique' in obj:
            obj['dejavu'] = set()

    # Parse for features of interest
    log.info('Reading %s...', input_file)
    with gzip.open(input_file) as input_fh:
        # Uncertain number of header lines, tricky to use csv.reader
        for line in input_fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            assert len(parts) == len(gtf_fields)
            feature_dict = dict(zip(gtf_fields, parts))
            if (feature_dict['feature'] in interest and
                feature_dict['seqname'] in seqname_dict):
                _handle(template, feature_dict, seqname_dict)

    # Output
    for key, obj in template.items():
        obj['output_fh'].close()
        log.info('Wrote %d items to %s', len(obj['dejavu']), obj['output_file'])


def _handle(template, feature_dict, seqname_dict):
    """
    Handle a single feature, write to file
    """

    # Deal with coordinates
    _transform_coordinates(feature_dict, seqname_dict)

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

            tup = tuple(map(model.get, obj['unique']))
            if tup in obj['dejavu']:
                raise Exception("Error %s occurs more than once" % str(tup))
            else:
                obj['dejavu'].add(tup)
                obj['writer'].writerow(model)


def _hash_dict(input_dict):
    inst = md5()
    inst.update(str(input_dict))
    return inst.hexdigest()


def db_transform(template, db_dict, require_ccds, namespace):
    """
    Transform data from raw parse to produce CSV files for database load
    """

    # Read gene file
    gene_dict = {}
    with open(template['gene']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for row_dict in reader:
            gene_dict[row_dict['accession']] = row_dict
    log.debug('Loaded %d genes', len(gene_dict))

    # Add a namespace for gene if not provided
    for gene_accession, gene in gene_dict.items():
        if not 'namespace' in gene or not gene['namespace']:
            gene['namespace'] = namespace

    # Read transcripts file
    transcript_dict = {}
    with open(template['transcript']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for row_dict in reader:
            transcript_dict[row_dict['accession']] = row_dict
    log.debug('Loaded %d transcripts', len(transcript_dict))

    # Read the coding exons and index conveniently on transcript_accession, transcript_order
    cds_dict = collections.OrderedDict()
    with open(template['CDS']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for exon_dict in reader:
            transcript_accession = exon_dict['transcript_accession']
            transcript_order = int(exon_dict['transcript_order'])
            if not transcript_accession in cds_dict:
                cds_dict[transcript_accession] = collections.OrderedDict()
            cds_dict[transcript_accession][transcript_order] = exon_dict
        input_fh.close()
    log.debug('Loaded %d CDS', len(cds_dict))

    # Read stop codons if present
    stop_dict = {}
    if 'stop_codon' in template:
        with open(template['stop_codon']['output_file'], 'r') as input_fh:
            reader = csv.DictReader(input_fh)
            for codon_dict in reader:
                transcript_accession = codon_dict['transcript_accession']
                transcript_order = int(codon_dict['transcript_order'])
                tup = transcript_accession, transcript_order
                if tup in stop_dict:
                    raise Exception('Multiple stop features for same codon %s', str(tup))
                stop_dict[tup] = codon_dict
            input_fh.close()
        log.debug('Loaded %d stop codons', len(stop_dict))

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

    # Calculate cds and cdsregions data
    #
    # The cds and cdsregion tables in the database represent
    # coding sequences i.e. the parts of a genomic sequence
    # that code for a protein, including the start and stop
    # codons. The cdsregion table holds coding exons. The
    # cds table is the transcript-level grouping. 
    #
    # Because transcripts have non-coding exons and exons 
    # which are partly coding and partly non-coding, multiple 
    # transcripts cans have the same coding sequence.
    #
    # In the following section we generate cds and cdsregion
    # records via a process of coordinate de-duplication.
    # We say that transcripts have the same coding sequence 
    # if their coding regions are on the same strand, in the 
    # same order and with the same coordinates. The code
    # identifies the distinct coding sequences and generates
    # the one to many mapping from coding sequence to 
    # transcript.
    #
    # The CCDS project 
    # http://www.ensembl.org/info/genome/genebuild/ccds.html
    # provides a curated set of coding sequences, for human
    # and mouse. Each coding sequence in their set carries
    # a CCDS id and the coordinates of the coding regions.
    # Further, the Ensembl GTF files for human and mouse 
    # include the CCDS id as an attribute of the transcript 
    # feature. 
    # 
    # The naive expectation would be that in those GTF files, 
    # transcripts with the same CCDS id should have the same
    # coding sequence, according to our definition above. 
    # In fact, this is almost always the case. The handful of
    # exceptions are logged. It should also be the case that 
    # equivalent transcripts should have the same CCDS id. 
    # Again this is almost always so. Exceptions are logged. 
    # Finally, one also might expect that every protein
    # coding transcript in the human and mouse GTF files 
    # should have a CCDS id, but this is not the case. There
    # are other 'quality criteria' used for the assignment of 
    # CCDS ids and about 1000 genes have no transcript with
    # a CCDS id.
    # 
    # Transcripts that have a coding sequence but no CCDS id
    # (e.g. not human or mouse, or simply no CCDS id assigned) 
    # are handled in this script according to the '--require-ccds' 
    # flag, as follows.
    #
    # If the --require-ccds flag is set true (default) the
    # coding sequence of a transcript with no CCDS id is 
    # omitted entirely and not represented in the cds and 
    # cdsregion files. 

    # If the --require-ccds flag is set false, all coding 
    # sequences are written out to the cds and cdsregion 
    # files.
    # 
    # Where CCDS id is available for a coding sequence, the
    # CCDS is is assigned to the 'name' field and the
    # 'is_consensus' flag is set to true. If there are 
    # alternative CCDS ids for the same coding sequence 
    # (unusual) separate cds record are generated for each 
    # instance. Thus the the cds records are not necessarily 
    # unique in terms of coordinates, but the duplication 
    # level is low and probably insignificant.
    #
    # If no CCDS id is available for a coding sequence, the
    # code attempts to used the 'protein_id' field from the
    # corresponding GTF records. If there are multiple 
    # protein ids, separate records are generated for each
    # instance. If there are no protein ids available for 
    # coding sequence, an auto-generated UID is assigned. 
    #
    # In all cases where a transcript has a coding sequence
    # but no CCDS id, the id assigned to the coding sequence
    # is copied back into the 'ccds_accession' field, to 
    # support front-end display of this information.

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
        # Possible non-uniqueness at coordinate level, but 
        # unusual for CCDS
        dict_to_hash = collections.OrderedDict({
            'name': None,
            'namespace': None,
            'is_consensus': 0,
            'exon_coords': []
        })
        for transcript_order in sorted(cds.keys()):
            exon_dict = cds[transcript_order]
            coord_dict = collections.OrderedDict({key: exon_dict[key] for key in (
                'chromosome_name', 'start', 'end', 'strand',
                'phase_start', 'phase_end',
                'cds_order')})
            dict_to_hash['exon_coords'].append(coord_dict)
        if not dict_to_hash['name'] and transcript_accession in transcript_to_ccds:
            dict_to_hash['name'] = transcript_to_ccds[transcript_accession]
            dict_to_hash['namespace'] = 'ccds'
            dict_to_hash['is_consensus'] = 1
        if not dict_to_hash['name'] and transcript_accession in cds_dict:
            cds = cds_dict[transcript_accession]
            for transcript_order, exon_dict in cds.items():
                if 'protein_id' in exon_dict:
                    dict_to_hash['name'] = exon_dict['protein_id']
                    break
        if not dict_to_hash['name']:
            dict_to_hash['name'] = 'CDS%06d' % (len(cds_to_transcript) + 1)
        hash_val = _hash_dict(dict_to_hash)
        if not hash_val in cds_to_transcript:
            cds_to_transcript[hash_val] = {
                'accession': "eq_{0}".format(len(cds_to_transcript)), 
                'namespace': dict_to_hash['namespace'],
                'name': dict_to_hash['name'],
                'gene_accession': transcript_dict[transcript_accession]['gene_accession'],
                'is_consensus': dict_to_hash['is_consensus'],
                'status': None,
                'transcript_accession_list': [],
                'exon_coords': dict_to_hash['exon_coords']
            }
            count += len(dict_to_hash['exon_coords'])
        cds_to_transcript[hash_val]['transcript_accession_list'].append(transcript_accession)
        transcript_to_cds[transcript_accession] = hash_val
    log.debug('Found %d coding sequences with total %d exons', 
              len(cds_to_transcript), count)

    # Calculate cds coordinates as the envelope of constituent exons
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

    # Check that transcripts with same CCDS have the same coding sequence
    for ccds_accession, transcript_accession_list in ccds_to_transcript.items():
        first_transcript_accession = transcript_accession_list[0]
        first_hash_val = transcript_to_cds[first_transcript_accession]
        for transcript_accession in transcript_accession_list:
            if not transcript_to_cds[transcript_accession] == first_hash_val:
                log.warning('Transcripts %s and %s with same CCDS have different coding sequences',
                            first_transcript_accession,
                            transcript_accession)

    # Assign namespace to transcripts
    for transcript_accession, transcript in transcript_dict.items():
        transcript['namespace'] = namespace

    # Write the cds name into the ccds_accession field of transcript,
    # if that field is not set - field may be displayed in front end
    for transcript_accession, transcript in transcript_dict.items():
        if not 'ccds_accession' in transcript or not transcript['ccds_accession']:
            if transcript_accession in transcript_to_cds:
                transcript['ccds_accession'] = cds_to_transcript[transcript_to_cds[transcript_accession]]['name']

    # Label each transcript with it's cds accession if available
    # This supports the fkey cds_id in transcript so you can tell
    # which transcripts belong to each equivalent set
    for transcript_accession, transcript in transcript_dict.items():
        transcript['cds_accession'] = None
        if transcript_accession in transcript_to_cds:
            transcript['cds_accession'] = cds_to_transcript[transcript_to_cds[transcript_accession]]['accession']

    # =============== OUTPUT ================== 

    # Produce the gene file
    fieldnames = [
        'accession',
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
        log.debug('Wrote %d regions to %s', count, db_dict['gene'])

    # Produce the cds file
    fieldnames = [
        'accession',
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
        log.debug('Wrote %d regions to %s', count, db_dict['cds'])

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
                coord_dict['cds_accession'] = obj['accession']
                row_dict = {key: coord_dict[key] for key in fieldnames}
                if int(row_dict['start']) > int(row_dict['end']):
                    log.warn('Omitting invalid cdsregion coordinates for CDS %s', obj['accession'])
                    continue
                writer.writerow(row_dict)
                count += 1
        output_fh.close()
        log.debug('Wrote %d exons to %s', count, db_dict['cdsregion'])

    # Produce the transcripts file
    fieldnames = [
        'accession',
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
        log.debug('Wrote %d regions to %s', count, db_dict['transcript'])

    # Read and write exon file - quite large so modify on the 
    # fly rather than holding in memory.
    # Merge the data for coding exons. Include the coding start and
    # end as separate columns, to support a possible future
    # coding_track_id in the exon table.
    # REGRESSION NOTE
    # In the current database, if the transcript has a cds entry,
    # the exon phase start/ends are OVERWRITTEN with the coding
    # coordinate start/end. The column 'phase_db_overwrite'
    # indicates to the loader when to do this.
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
        with open(template['exon']['output_file'], 'r') as input_fh:
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
        log.debug('Wrote %d exons to %s', count, db_dict['exon'])


def _transform_coordinates(feature_dict, seqname_dict):

    # Sequence name must be in the config
    if not feature_dict['seqname'] in seqname_dict:
        raise Exception('Invalid seqname, does not match config')
    feature_dict['seqname'] = seqname_dict[feature_dict['seqname']]

    # Strand must be '+' or '-'
    if not feature_dict['strand'] in ('+', '-'):
        raise Exception('Invalid strand')
    feature_dict['strand'] = 1 if feature_dict['strand'] == '+' else -1

    # Input is one-based [,] per GTF specification
    # CSV output is zero-based [,)
    feature_dict['start'] = str(int(feature_dict['start']) - 1)


if __name__ == '__main__':
    main()

