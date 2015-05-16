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
# database annotations in the database. 
#
# Logically the second stage could be in a separate script,
# but practically the two stages will almost always be run
# in sequence, so the two are together. The second stage 
# depends on particular values for the aliases in the JSON 
# template.

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
    cdsregion_fieldnames = None
    with open(template['CDS']['output_file'], 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for exon_dict in reader:
            if not cdsregion_fieldnames:
                cdsregion_fieldnames = reader.fieldnames
            transcript_accession = exon_dict['transcript_accession']
            transcript_order = exon_dict['transcript_order']
            if not transcript_accession in cds_dict:
                cds_dict[transcript_accession] = collections.OrderedDict()
            cds_dict[transcript_accession][int(transcript_order)] = exon_dict
        input_fh.close()
    log.debug('Loaded %d CDS', len(cds_dict))

    # Read stop codons if present
    stop_dict = {}
    if 'stop_codon' in template:
        with open(template['stop_codon']['output_file'], 'r') as input_fh:
            reader = csv.DictReader(input_fh)
            for codon_dict in reader:
                tup = codon_dict['transcript_accession'], int(codon_dict['transcript_order'])
                # Already checked for uniqueness
                stop_dict[tup] = codon_dict
            input_fh.close()
        log.debug('Loaded %d stop codons', len(stop_dict))

    if len(stop_dict):
        # Append stop codons
        # Note there are sometimes split codongs (length % 3 ! = 0) so 
        # important to do this before calculating phases
        for transcript_accession, cds in cds_dict.items():
            for transcript_order, exon_dict in cds.items():
                tup = transcript_accession, transcript_order
                # Append stop codon
                if tup in stop_dict:
                    codon_dict = stop_dict[tup]
                    exon_dict['start'] = min(int(exon_dict['start']),
                                             int(codon_dict['start']))
                    exon_dict['end'] = max(int(exon_dict['end']),
                                           int(codon_dict['end']))

    # Calculate phases
    # Frame is the offset of the 1st codon from the feature start:
    # http://www.ensembl.org/info/website/upload/gff.html
    # Phase is the offset of the feature start from the 1st codon:
    # https://github.com/Ensembl/ensembl/blob/release/79/modules/Bio/EnsEMBL/Exon.pm
    # Hence phase = ( - frame ) % 3
    for transcript_accession, cds in cds_dict.items():
        for transcript_order, exon_dict in cds.items():
            exon_length = int(exon_dict['end']) - int(exon_dict['start'])
            exon_dict['phase_start'] = (- int(exon_dict['frame'])) % 3
            exon_dict['phase_end'] = (exon_dict['phase_start'] + exon_length) % 3

    # Assign cds_order
    # We could use the transcript_order for ordering, but that
    # doesn't always start from 1
    for transcript_accession, cds in cds_dict.items():
        order = 1
        for transcript_order in sorted(cds.keys()):
            cds[transcript_order]['cds_order'] = order
            order += 1

    # Calculate cds and cdsregions data
    #
    # The cds and cdsregion tables in the database hold
    # transcript and exon level structures that are derived
    # from GTF transcript and exon records through a process
    # of coordinate de-duplication.
    #
    # Two transcripts are equivalent (duplicates) if they are
    # are on the same strand, have the same number of exons
    # and if the coordinates of the exons are identical.
    #
    # The code calculates the sets of equivalent transcripts.
    # The output cds file lists one record for each equivalent 
    # set. The cdsregion file holds one set of coordinates 
    # for the corresponding constituent exons.
    #
    # For human and mouse genomes, the CCDS project 
    # (http://www.ensembl.org/info/genome/genebuild/ccds.html)
    # provides just such a non-redundant set of transcripts.
    # Equivalent transcrpts (according to CCDS) are labelled
    # with the same CCDS id. Further, the Ensembl GTF files
    # for human and mouse include the CCDS id as an attribute
    # of the transcript feature. So the expectation would be
    # that in those files, transcripts with the same CCDS id 
    # should be equivalent (according to our definition). This 
    # is almost always the case. Exceptions are logged. It
    # should also be the case that equivalent transcripts
    # should have the same CCDS id. Again this is almost always
    # so. Exceptions are logged. Finally, one also might 
    # expect that every equivalent set in the human and
    # mouse GTF files should have a CCDS id, but this is not
    # the case because there are other 'quality criteria' used
    # for the assignment of CCDS ids. 
    # 
    # Equivalent sets with no CCDS id (e.g. not human or mouse,
    # or simply no CCDS id assigned) are handled according to 
    # a '--require-ccds' flag, as follows.
    #
    # If the --require-ccds flag is set true (default) equivalent
    # sets with with no CCDS id are ignored, omitted entirely
    # from the cds and cdsregion files. In this case the name
    # field of the cds file holds the CCDS id, labelling the
    # equivalent set and the 'is_consensus' field is set True.
    #
    # If the --require-ccds flag is set false, all equivalent
    # sets are represented in the cds and cdsregion files.
    # The name field in the cds file uses the CCDS id if it is
    # available. Otherwise it attempts to use the 'protein_id'
    # field from a corresponding GTF record. If this is
    # inconsistent or null, it assigns an auto-generated UID.
    # All assigned names are written out in the 'ccds_accession'
    # field of the transcript file. So in this case the
    # 'ccds_accession' field represents our own definition of
    # a 'consensus CDS' rather than an official CCDS id.

    # Maps to enable consistency check for CCDS, where provided
    ccds_to_transcript = {}
    transcript_to_ccds = {}
    for transcript_accession, transcript in transcript_dict.items():
        if 'ccds_accession' in transcript and transcript['ccds_accession']:
            ccds_accession = transcript['ccds_accession']
            if not ccds_accession in ccds_to_transcript:
                ccds_to_transcript[ccds_accession] = []
            ccds_to_transcript[ccds_accession].append(transcript_accession)
            transcript_to_ccds[transcript_accession] = ccds_accession

    # Calculate equivalent sets
    cds_to_transcript = {}
    transcript_to_cds = {}
    cdsregion_count = 0
    for transcript_accession, cds in cds_dict.items():
        # Ignore transcripts without a ccds_id, if the require_ccds flag is set
        if require_ccds:
            if not transcript_accession in transcript_to_ccds:
                continue
        exon_coords = []
        for transcript_order in sorted(cds.keys()):
            exon_dict = cds[transcript_order]
            coord_dict = {key: exon_dict[key] for key in (
                'chromosome_name', 'start', 'end', 'strand',
                'phase_start', 'phase_end',
                'cds_order')}
            exon_coords.append(coord_dict)
        hash_val = _hash_dict({'exon_coords': exon_coords})
        if not hash_val in cds_to_transcript:
            cds_to_transcript[hash_val] = {
                'accession': "eq_{0}".format(len(cds_to_transcript)),
                'namespace': None,
                'name': None,
                'gene_accession': transcript_dict[transcript_accession]['gene_accession'],
                'is_consensus': 0 if require_ccds else 1,
                'status': None,
                'transcript_accession_list': [],
                'exon_coords': exon_coords
            }
            cdsregion_count += len(exon_coords)
        cds_to_transcript[hash_val]['transcript_accession_list'].append(transcript_accession)
        transcript_to_cds[transcript_accession] = hash_val
    log.debug('Found %d sets of equivalent transcripts (cds)', len(cds_to_transcript))
    log.debug('Corresponding to %d exons (cdsregions)', cdsregion_count)

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
        obj['start'] = start
        obj['end'] = end

    # Check that transcripts with same CCDS are equivalent
    for ccds_accession, transcript_accession_list in ccds_to_transcript.items():
        first_transcript_accession = transcript_accession_list[0]
        first_hash_val = transcript_to_cds[first_transcript_accession]
        for transcript_accession in transcript_accession_list:
            if not transcript_to_cds[transcript_accession] == first_hash_val:
                log.warning('Transcripts %s and %s with same CCDS are not equivalent',
                            first_transcript_accession,
                            transcript_accession)

    # Check that equivalent transcripts have the same CCDS, if they have one
    for hash_val, obj in cds_to_transcript.items():
        first_transcript_accession = obj['transcript_accession_list'][0]
        if 'ccds_accession' in transcript_to_ccds:
            first_ccds_accession = transcript_to_ccds[first_transcript_accession]
            for transcript_accession in obj['transcript_accession_list']:
                if 'ccds_accession' in transcript_to_ccds:
                    if not transcript_to_ccds[transcript_accession] == first_ccds_accession:
                        log.warn('Equivalent transcripts %s and %s have different CCDS ids',
                                 first_transcript_accession,
                                 transcript_accession)

    # Assign the CCDS id as the cds name if we have one
    for hash_val, obj in cds_to_transcript.items():
        for transcript_accession in obj['transcript_accession_list']:
            if transcript_accession in transcript_to_ccds:
                obj['name'] = transcript_to_ccds[transcript_accession]
                obj['is_consensus'] = 1
                obj['namespace'] = 'ccds'
                obj['status'] = 'Public'

    # Try to assign protein id as a cds name
    for hash_val, obj in cds_to_transcript.items():
        first_transcript_accession = None
        first_protein_accession = None
        consistent = True
        for transcript_accession in obj['transcript_accession_list']:
            if not consistent:
                break
            cds = cds_dict[transcript_accession]
            for transcript_order, exon_dict in cds.items():
                if not consistent:
                    break
                if 'protein_id' in exon_dict:
                    if not first_protein_accession:
                        first_transcript_accession = transcript_accession
                        first_protein_accession = exon_dict['protein_id']
                    if not first_protein_accession == exon_dict['protein_id']:
                        consistent = False
                        log.warn('Inconsistent protein ids for transcripts %s %s',
                                 first_transcript_accession,
                                 transcript_accession)
        if consistent and first_protein_accession:
            obj['name'] = first_protein_accession

    # Assign a UID for any un-named remainders
    uid_counter = 1
    for hash_val, obj in cds_to_transcript.items():
        if not obj['name']:
            obj['name'] = 'EQS%06d' % (uid_counter)
            uid_counter += 1

    # Write the cds name into the ccds_accession field of transcript,
    # if it is not already there
    for transcript_accession, transcript in transcript_dict.items():
        transcript['namespace'] = namespace
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
                        # This is a coding exon
                        coding_exon_dict = cds[transcript_order]
                        exon_dict['phase_start'] = coding_exon_dict['phase_start']
                        exon_dict['phase_end'] = coding_exon_dict['phase_end']
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

