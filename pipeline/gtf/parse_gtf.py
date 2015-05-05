"""Extract Gene, Transcript, Exon, CDS features from EMBL GTF format file"""
import gzip
import re
import csv
import gzip
import logging
import click
import yaml
import collections

from pprint import pprint

logging.basicConfig()
log = logging.getLogger(__name__)    
log.setLevel(logging.INFO)

@click.group()
def cli():
    pass

def _base_dict(id_key, id_dict, coord_dict):

    base = collections.OrderedDict()

    if id_key and id_dict:
        base[id_key] = id_dict[id_key]

    if coord_dict:
        # Conform to genomebrowser convention
        base['chromosome_name'] = coord_dict['seqname']
        base['strand'] = ( 1 if coord_dict['strand'] == '+' else -1 )
        base['start'] = int(coord_dict['start'])
        base['end'] = int(coord_dict['end'])

    return base


@cli.command()
@click.argument('embl_gtf')
def parse_gtf(embl_gtf):
    """
    Extract genes/transcripts/exons from EMBL GTF
    """

    # Reference
    # http://www.ensembl.org/info/website/upload/gff.html

    gtf_fields = ['seqname',
                  'source',
                  'feature',
                  'start',
                  'end',
                  'score',
                  'strand',
                  'frame',
                  'attribute']

    features_of_interest = ['transcript',
                            'exon',
                            'CDS']
    
    # Parse for features of interest
    gtf_list = []

    log.info('Reading %s...', embl_gtf)
    with gzip.open(embl_gtf) as input_fh:
        for line in input_fh:
            if not line.startswith('#'):
                break
        csv_reader = csv.reader(input_fh, delimiter='\t')
        for row in csv_reader:
            gtf_dict = collections.OrderedDict(zip(gtf_fields, row))
            if gtf_dict['feature'] in features_of_interest:
                gtf_list.append(gtf_dict)
#                if len(gtf_list) >= 10000:
#                    break

    log.info('Extracted %d features', len(gtf_list)) 

    # Unpack the attribute field of each feature into a second level dict
    attribute_regex = r"(.*?)\s+\"(.*?)\";\s+"                               
    for gtf_dict in gtf_list:
        if gtf_dict['attribute']:
            pairs = re.findall(attribute_regex, gtf_dict['attribute'])
            # Overwrite the string with the dict!
            gtf_dict['attribute'] = collections.OrderedDict(pairs)

    # Collect genes
    gene_dict = collections.OrderedDict()
    for gtf_dict in gtf_list:
        if gtf_dict['feature'] == 'transcript':
            # Base id, no coordinates
            gene = _base_dict('gene_id', gtf_dict['attribute'], None)
            # Copy gene related keys
            for key in gtf_dict['attribute']:
                if key.startswith('gene_') and key not in gene:
                    gene[key] = gtf_dict['attribute'][key]
            # Append to list unless already there
            if not gene['gene_id'] in gene_dict:
                gene_dict[gene['gene_id']] = gene
            else:
                # Check identity
                assert not cmp(gene_dict[gene['gene_id']], gene)

    # Collect transcripts, label by gene
    transcript_dict = collections.OrderedDict()
    for gtf_dict in gtf_list:
        if gtf_dict['feature'] == 'transcript':
            # Base id and coordinates
            transcript = _base_dict('transcript_id', gtf_dict['attribute'], gtf_dict)
            # Copy transcript related keys
            for key in gtf_dict['attribute']:
                if key.startswith('transcript_') and key not in transcript:
                    transcript[key] = gtf_dict['attribute'][key]
            # Insert gene ID
            gene_id = gtf_dict['attribute']['gene_id']
            if not gene_id in gene_dict:
                raise Exeception('Could not find %s', gene_id)
            transcript['gene_id'] = gene_id
            # List for exons (see below)
            transcript['exons'] = []
            # Append to list unless already there
            if not transcript['transcript_id'] in transcript_dict:
                transcript_dict[transcript['transcript_id']] = transcript
            else:
                raise Exception('Already have %s', transcript['transcript_id'])

    # An exon feature in the GTF file carries an exon_id, 
    # exon_number, transcript_id and in some cases an 
    # exon_version. The same exon_id can appear with 
    # different transcript_ids, and in general when it 
    # does so, the exon_number and exon_versio are 
    # specific to that occurrence. 
    # The code asserts that exon features with the same 
    # exon_id must have the same coordinates.
    # To represent the relationship in a non-redundant
    # way, this code represents an exon with exon_id and
    # coordinates, and the transcript dict holds a list 
    # exon-tuples representing the occcurence of specified
    # exons in the transcript. An exon-tuple contains 
    # exon_id, exon_number (transcript specific) and
    # (in some cases) a third item the exon_version 
    # (assumed to be transcript specific). 

    # Collect exons, append to transcripts
    exon_dict = {}
    for gtf_dict in gtf_list:
        if gtf_dict['feature'] == 'exon':
            # Base id and coordinates
            exon = _base_dict('exon_id', gtf_dict['attribute'], gtf_dict)
            # Copy all exon related keys (except exon_number)
            for key in gtf_dict['attribute']:
                # Ignore some keys
                if key in ['exon_number', 'exon_version']:
                    continue
                if key.startswith('exon_') and key not in exon:
                    exon[key] = gtf_dict['attribute'][key]
            # Append exon_id, exon_number tuple to transcript
            transcript_id = gtf_dict['attribute']['transcript_id'] 
            if transcript_id not in transcript_dict:
                raise Exception('Could not find %s' % (transcript_id))
            # Tuple
            if 'exon_version' in gtf_dict['attribute']:
                transcript_dict[transcript_id]['exons'].append(
                    (exon['exon_id'], gtf_dict['attribute']['exon_number'], gtf_dict['attribute']['exon_version']))
            else:
                transcript_dict[transcript_id]['exons'].append(
                    (exon['exon_id'], gtf_dict['attribute']['exon_number']))
            # Append to list unless already there
            if not exon['exon_id'] in exon_dict:
                exon_dict[exon['exon_id']] = exon
            else:
                # Check identity 
                assert not cmp(exon_dict[exon['exon_id']], exon)

    # Collect cdsregions, label by transcript
    cdsregion_dict = collections.OrderedDict()
    for gtf_dict in gtf_list:
        if gtf_dict['feature'] == 'CDS':
            # No id, coordinates
            cdsregion = _base_dict(None, None, gtf_dict)
            # Collect
            for key in 'protein_id', 'exon_number':
                cdsregion[key] = gtf_dict['attribute'][key];
            # Insert transcript ID
            transcript_id = gtf_dict['attribute']['transcript_id'] 
            if transcript_id not in transcript_dict:
                raise Exception('Could not find %s' % (transcript_id))
            cdsregion['transcript_id'] = transcript_id
            # Append to list unless already there
            check_for = cdsregion['protein_id'], cdsregion['exon_number']
            if not check_for in cdsregion_dict:
                cdsregion_dict[check_for] = cdsregion
            else:
                raise Exception('Already have %s' % (" ".join(check_for)))

    # The exon, gene and CDS dicts are all flat, so 
    # output as CSV.
    for output_file, obj_dict in [
            ('gene', gene_dict),
            ('exon', exon_dict),
            ('cdsregion', cdsregion_dict)]:
        with open(output_file + '.csv', 'w') as output_fh:
            done_header = 0
            writer = None
            for item_id in obj_dict:
                if not writer:
                    fieldnames = obj_dict[item_id].keys()
                    writer = csv.DictWriter(output_fh, fieldnames)
                    writer.writeheader()
                writer.writerow(obj_dict[item_id])
            output_fh.close()
            log.info('Wrote %d items to %s', len(obj_dict), output_file)

    # The transcript dict contains a list of exon-tuples
    # so output in YAML.
    for output_file, obj_dict in [
            ('transcript', transcript_dict)]:
        with open(output_file + '.yaml', 'w') as output_fh:
            output_fh.write(yaml.dump(obj_dict))
            output_fh.close()
            log.info('Wrote %d items to %s', len(obj_dict), output_file)

if __name__ == '__main__':
    cli()

