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

    print("COMPARING %s -------" % (title))
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
    checked = 0
    for key in gtf:
        if key in db:
            gtf_dict =reg_dict(gtf[key])
            db_dict = reg_dict(db[key])
            if cmp(gtf_dict, db_dict):
                print('Difference!')
                print('GTF: ' + str(gtf_dict))
                print('DB:  ' + str(db_dict))
                count += 1
    print("Total %d checked" % (checked))
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
            'namespace': gene.namespace,
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
            'namespace': transcript.namespace,
            'biotype': transcript.biotype,
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
            # Fields not in db
            transcript_dict.pop('transcript_support_level')
            transcript_dict.pop('cds_accession')
            if transcript_dict['chromosome_name'] == chromosome.name:
                gtf_transcripts_dict[transcript_dict['accession']] = transcript_dict
        
    log.info('Read %d transcripts from %s', len(gtf_transcripts_dict), input_file)

    compare_lists(db_transcripts_dict, gtf_transcripts_dict, 'transcripts')

def compare_cdsregions(chromosome, cds_map, input_file):
    """ 
    """     

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    # Fetch the cdsregions from the database for the given chromosome
    db_cdsregions_dict = {}
    for cdsregion in (S.query(gbm.CdsRegion)
                      .join(gbm.Coordinates)
                      .join(gbm.Chromosome)
                      .filter(gbm.Chromosome.id == chromosome.id)
                      .all()):
        cdsregion_dict = {
            'cds_accession': cds_map[cdsregion.cds_id],
            'chromosome_name': chromosome.name,
            'start': cdsregion.coordinates.start_end.lower,
            'end': cdsregion.coordinates.start_end.upper,
            'strand': cdsregion.coordinates.strand,
            'cds_order': cdsregion.cds_order,
            'phase_start': str(cdsregion.phase_start),
            'phase_end': str(cdsregion.phase_end)
        }
        cdsregion_dict = reg_dict(cdsregion_dict)
        tup = (cdsregion_dict['cds_accession'],
               int(cdsregion_dict['cds_order']))
        db_cdsregions_dict[tup] = cdsregion_dict
        
    log.info('Read %d cdsregions from database', len(db_cdsregions_dict))

    # Read the CSV file
    gtf_cdsregions_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for cdsregion_dict in reader:
            if cdsregion_dict['chromosome_name'] == chromosome.name:
                tup = (cdsregion_dict['cds_accession'],
                       int(cdsregion_dict['cds_order']))
                gtf_cdsregions_dict[tup] = cdsregion_dict
        
    log.info('Read %d cdsregions from %s', len(gtf_cdsregions_dict), input_file)

    compare_lists(db_cdsregions_dict, gtf_cdsregions_dict, 'cdsregions')


def compare_cds(chromosome, input_file):
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
    db_cds_dict = {}
    # Need this for later check on cdsregions
    cds_map = {}
    for cds in (S.query(gbm.CDS)
                .join(gbm.Coordinates)
                .join(gbm.Chromosome)
                .filter(gbm.Chromosome.id == chromosome.id)
                .all()):
        cds_dict = {
            'accession': cds.id,  # Overwritten later
            'name': cds.name,
            'namespace': cds.namespace,
            'is_consensus': int(cds.is_consensus),
            'gene_accession': genes_revmap[cds.gene_id],
            'chromosome_name': chromosome.name,
            'start': cds.coordinates.start_end.lower,
            'end': cds.coordinates.start_end.upper,
            'strand': cds.coordinates.strand,
            'status': cds.status
        }
        tup = (cds.name,
               str(cds.coordinates.start_end.lower),
               str(cds.coordinates.start_end.upper),
               str(cds.coordinates.strand))            
        if tup in db_cds_dict:
            raise Exception('Invalid tup, not unique')
        cds_map[tup] = cds.id
        db_cds_dict[tup] = cds_dict
        
    log.info('Read %d cds from database', len(db_cds_dict))

    # Read the CSV file
    gtf_cds_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for cds_dict in reader:
            if cds_dict['chromosome_name'] == chromosome.name:
                tup = (cds_dict['name'],
                       cds_dict['start'],
                       cds_dict['end'],
                       cds_dict['strand'])
                if tup in gtf_cds_dict:
                    raise Exception('Invalid tup, not unique')
                if tup in db_cds_dict:
                    db_cds_dict[tup]['accession'] = cds_dict['accession']
                gtf_cds_dict[tup] = cds_dict
                cds_map[cds_map[tup]] = cds_dict['accession']
        
    log.info('Read %d cds from %s', len(gtf_cds_dict), input_file)

    compare_lists(db_cds_dict, gtf_cds_dict, 'cds')
    
    return cds_map
    
def compare_exons(chromosome, input_file):
    """ 
    """     

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    transcript_revmap = {}
    for transcript in (S.query(gbm.Transcript)
                       .join(gbm.Coordinates)
                       .join(gbm.Chromosome)
                       .filter(gbm.Chromosome.id == chromosome.id)
                       .all()):
        transcript_revmap[transcript.id] = transcript.accession
        
    # Fetch the cdsregions from the database for the given chromosome
    db_exon_dict = {}
    for exon in (S.query(gbm.Exon)
                 .join(gbm.Coordinates)
                 .join(gbm.Chromosome)
                 .filter(gbm.Chromosome.id == chromosome.id)
                 .all()):
        exon_dict = {
            'transcript_accession': transcript_revmap[exon.transcript_id],
            'chromosome_name': chromosome.name,
            'start': exon.coordinates.start_end.lower,
            'end': exon.coordinates.start_end.upper,
            'strand': exon.coordinates.strand,
            'transcript_order': exon.transcript_order,
            'phase_start': str(exon.phase_start),
            'phase_end': str(exon.phase_end)
        }
        exon_dict = reg_dict(exon_dict)
        tup = (exon_dict['transcript_accession'],
               int(exon_dict['transcript_order']))
        db_exon_dict[tup] = exon_dict
        
    log.info('Read %d exons from database', len(db_exon_dict))

    # Read the CSV file
    gtf_exon_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for exon_dict in reader:
            if exon_dict['chromosome_name'] == chromosome.name:
                # Fields not in db
                exon_dict.pop('coding_start')
                exon_dict.pop('coding_end')
                exon_dict.pop('phase_db_overwrite')
                tup = (exon_dict['transcript_accession'],
                       int(exon_dict['transcript_order']))
                gtf_exon_dict[tup] = exon_dict
        
    log.info('Read %d exons from %s', len(gtf_exon_dict), input_file)

    compare_lists(db_exon_dict, gtf_exon_dict, 'exons')




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
        cds_id_to_accession = compare_cds(chromosome, cds_file)
        compare_cdsregions(chromosome, cds_id_to_accession, cdsregions_file)
    
if __name__ == "__main__":
    main()
