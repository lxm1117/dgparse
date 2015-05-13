"""Compare output from GTF parser with database

Usage:
  regression.py (--genome=<name>) (--genes=<file>) (--transcripts=<file>) (--exons=<file>) (--cds=<file>) (--cdsregions=<file>) [--log=<loglevel>]

Options:
  --help                    Show this screen
  --genome=<name>           Genome version (e.g. GRCh38.p2)
  --genes=<file>            CSV file of genes parsed from GTF
  --transcripts=<file>      CSV file of transcripts parsed from GTF
  --exons=<file>            CSV file of exons parsed from GTF
  --cds=<file>              CSV file of cds data parsed from GTF
  --cdsregions=<file>       CSV file of cdsregions parsed from GTF
  --log=<loglevel>          Logging level [default: ERROR]

"""
import os
import logging
import subprocess
import re
import csv
import collections

from docopt import docopt
from pprint import pprint

import addresses
from manifest.fulfilment.common import session
import genomebrowser.models as gbm

from bio.sequtils import get_reverse_complement

logging.basicConfig()
log = logging.getLogger(__name__)

def reg_dict(my_dict):
    kvp = collections.OrderedDict()
    for key in sorted(my_dict.keys()):
        if my_dict[key]:
            kvp[key] = str(my_dict[key])
        else:
            kvp[key] = ''
    return kvp


def compare_lists(db, gtf, title):

    print("%s missing from DB?" % (title))
    count = 0
    for key in gtf:
        if not key in db:
            pprint(gtf[key])
            count += 1
    print("Total %d" % (count))
    log.info('Total %d missing from DB', count)
    count = 0
    print("%s missing from GTF?" % (title))
    for key in db:
        if not key in gtf:
            pprint(db[key])
            count += 1
    print("Total %d" % (count))
    log.info('Total %d missing from GTF', count)

    print("Comparing common set...")
    count = 0
    for key in gtf:
        if key in db:
            gtf_dict =reg_dict(gtf[key])
            db_dict = reg_dict(db[key])
            if cmp(gtf_dict, db_dict):
                print('Difference!')
                print('GTF: ' + str(gtf_dict))
                print('DB:  ' + str(db_dict))
                count += 1
    print("Total %d differences" % (count))
    log.info("Total %d differences", count)


def compare_genes(chromosome, input_file):
    """ 
    """     

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    # Fetch the genes from the database for the given chromosome
    db_genes_dict = {}
    for gene in (S.query(gbm.Gene)
                 .join(gbm.Coordinates)
                 .join(gbm.Chromosome)
                 .filter(gbm.Chromosome.id == chromosome.id)
                 .all()):
        gene_dict = {
            'accession': gene.accession,
            'biotype': gene.biotype,
            'name': gene.name,
            'chromosome_name': chromosome.name,
            'start': gene.coordinates.start_end.lower,
            'end': gene.coordinates.start_end.upper,
            'strand': gene.coordinates.strand
        }
        db_genes_dict[gene_dict['accession']] = gene_dict
    
    log.info('Read %d genes from database', len(db_genes_dict))

    # Read the CSV file
    gtf_genes_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for gene_dict in reader:
            if gene_dict['chromosome_name'] == chromosome.name:
                gtf_genes_dict[gene_dict['accession']] = gene_dict
        
    log.info('Read %d genes from %s', len(gtf_genes_dict), input_file)

    compare_lists(db_genes_dict, gtf_genes_dict, 'genes')


def compare_transcripts(chromosome, input_file):
    """ 
    """     

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    genes_revmap = {}
    for gene in (S.query(gbm.Gene)
                 .join(gbm.Coordinates)
                 .join(gbm.Chromosome)
                 .filter(gbm.Chromosome.id == chromosome.id)
                 .all()):
        genes_revmap[gene.id] = gene.accession

    # Fetch the genes from the database for given chromosome
    db_transcripts_dict = {}
    for transcript in (S.query(gbm.Transcript)
                 .join(gbm.Coordinates)
                 .join(gbm.Chromosome)
                 .filter(gbm.Chromosome.id == chromosome.id)
                 .all()):
        transcript_dict = {
            'gene_accession': genes_revmap[transcript.gene_id],
            'name': transcript.name,
            'accession': transcript.accession,
            'ccds_accession': transcript.ccds_accession,
            'chromosome_name': chromosome.name,
            'start': transcript.coordinates.start_end.lower,
            'end': transcript.coordinates.start_end.upper,
            'strand': transcript.coordinates.strand
        }
        db_transcripts_dict[transcript_dict['accession']] = transcript_dict
    
    log.info('Read %d transcripts from database', len(db_transcripts_dict))

    # Read the CSV file
    gtf_transcripts_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for transcript_dict in reader:
            if transcript_dict['chromosome_name'] == chromosome.name:
                gtf_transcripts_dict[transcript_dict['accession']] = transcript_dict
        
    log.info('Read %d transcripts from %s', len(gtf_transcripts_dict), input_file)

    compare_lists(db_transcripts_dict, gtf_transcripts_dict, 'transcripts')


def compare_cdsregions(chromosome, input_file):
    """ 
    """     

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    transcripts_revmap = {}
    ccds_map = {}
    for transcript in (S.query(gbm.Transcript)
                       .join(gbm.Coordinates)
                       .join(gbm.Chromosome)
                       .filter(gbm.Chromosome.id == chromosome.id)
                       .all()):
        transcripts_revmap[transcript.id] = transcript.accession
        if transcript.ccds_accession:
            ccds_map[transcript.ccds_accession] = transcript.accession

    cds_revmap = {}
    for cds in (S.query(gbm.CDS)
                .join(gbm.Coordinates)
                .join(gbm.Chromosome)
                .filter(gbm.Chromosome.id == chromosome.id)
                .all()):
        cds_revmap[cds.id] = cds.name

    # Fetch the cdsregions from the database for the given chromosome
    db_cdsregions_dict = {}
    for cdsregion in (S.query(gbm.CdsRegion)
                      .join(gbm.Coordinates)
                      .join(gbm.Chromosome)
                      .filter(gbm.Chromosome.id == chromosome.id)
                      .all()):
        ccds_accession = cds_revmap[cdsregion.cds_id]
        cdsregion_dict = {
            'transcript_name': ccds_map[ccds_accession],
            # Suspect NCBI has cds_order the other way around
            'chromosome_name': chromosome.name,
            'start': cdsregion.coordinates.start_end.lower,
            'end': cdsregion.coordinates.start_end.upper,
            'strand': cdsregion.coordinates.strand,
            'cds_order': cdsregion.cds_order,
            'phase_start': str(cdsregion.phase_start),
            #'phase_end': str(cdsregion.phase_end)
        }
        cdsregion_dict = reg_dict(cdsregion_dict)
        db_cdsregions_dict[(cdsregion_dict['transcript_name'],
                            cdsregion_dict['start'])] = cdsregion_dict
        
    log.info('Read %d cdsregions from database', len(db_cdsregions_dict))

    # Read the CSV file
    gtf_cdsregions_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for cdsregion_dict in reader:
            if cdsregion_dict['chromosome_name'] == chromosome.name:
                cdsregion_dict.pop('transcript_order')
                cdsregion_dict.pop('phase_end')
                cdsregion_dict.pop('frame')
                gtf_cdsregions_dict[(cdsregion_dict['transcript_name'],
                                     cdsregion_dict['start'])] = cdsregion_dict
        
    log.info('Read %d cdsregions from %s', len(gtf_cdsregions_dict), input_file)

    compare_lists(db_cdsregions_dict, gtf_cdsregions_dict, 'cdsregions')


def main():

    # Parse args
    args = docopt(__doc__, version='Breakpoints Annotator')

    # Set logging level
    loglevel = getattr(logging, args['--log'].upper())
    log.setLevel(loglevel)

    # Log args
    log.info(str(args))

    # Genome version
    genome_version = args['--genome']

    # Output file
    genes_file = os.path.abspath(args['--genes'])
    transcripts_file = os.path.abspath(args['--transcripts'])
    exons_file = os.path.abspath(args['--exons'])
    cds_file = os.path.abspath(args['--cds'])
    cdsregions_file = os.path.abspath(args['--cdsregions'])
    
    # Split processing by chromosome - large amount of data
    S = session(addresses.genomebrowser.database.name)
    chromosome_list = (S.query(gbm.Chromosome)
                       .join(gbm.Genome)
                       .filter_by(version=genome_version)
                       .order_by(gbm.Chromosome.length)
                       .all())
    
    for chromosome in chromosome_list:
        if chromosome.name != 'chr21':
            print('Warning restricting test to specific chromosome')
            continue
        log.info('===== Comparisons for %s =====', chromosome.name)
        compare_genes(chromosome, genes_file)
        compare_transcripts(chromosome, transcripts_file)
        compare_exons(chromosome, exons_file)
        compare_cds(chromosome, cds_file)
        compare_cdsregions(chromosome, cdsregions_file)
 
if __name__ == "__main__":
    main()
