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
import json

from operator import itemgetter
from docopt import docopt
from hashlib import md5
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
        if not my_dict[key] is None:
            kvp[key] = str(my_dict[key])
        else:
            kvp[key] = ''
    return kvp


def _hash_dict(input_dict):
#    return(str(input_dict))
    inst = md5()
    inst.update(str(input_dict))
    return inst.hexdigest()


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

    gtf_tmp = "{0}-{1}.tmp".format(title, 'gtf')
    with open(gtf_tmp, 'w') as tmp_fh:
        for key in sorted(gtf.keys()):
            tmp_fh.write("%s, %s\n" % (key, str(reg_dict(gtf[key]))))
    log.info('Wrote gtf data to %s', gtf_tmp)
    db_tmp = "{0}-{1}.tmp".format(title, 'db')
    with open(db_tmp, 'w') as tmp_fh:
        for key in sorted(db.keys()):
            tmp_fh.write("%s, %s\n" % (key, str(reg_dict(db[key]))))
    log.info('Wrote db data to %s', db_tmp)

    print("Comparing common set...")
    count = 0
    checked = 0
    for key in gtf:
        if key in db:
            gtf_dict = reg_dict(gtf[key])
            db_dict = reg_dict(db[key])
            if cmp(gtf_dict, db_dict):
                print '*************** %s difference...' % (title)
                print 'GTF: ' + str(gtf_dict)
                print 'DB:  ' + str(db_dict)
                count += 1
    print("Total %d checked" % (checked))
    msg = "Total %d %s differences" % (count, title)
    if count == 0:
        msg += ' *** PASS'
    else:
        msg += ' *** CHECK'
    print msg
    log.info(msg)


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

    compare_lists(db_genes_dict, gtf_genes_dict, 'gene')


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
            # Extra ';' in 45 transcripts, disable for now
            #'name': transcript.name,
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
            # Extra ';' in 45 transcripts, disable for now
            transcript_dict.pop('name')
            if transcript_dict['chromosome_name'] == chromosome.name:
                gtf_transcripts_dict[transcript_dict['accession']] = transcript_dict
        
    log.info('Read %d transcripts from %s', len(gtf_transcripts_dict), input_file)

    compare_lists(db_transcripts_dict, gtf_transcripts_dict, 'transcript')

def create_hash_to_cds_id_map(chromosome):
    """ 
    """     

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    db_cds_dict = {}
    for cds in (S.query(gbm.CDS)
                .join(gbm.Coordinates)
                .join(gbm.Chromosome)
                .filter(gbm.Chromosome.id == chromosome.id)
                .all()):
        cds_dict = {
            'name': cds.name,
            'namespace': cds.namespace,
            'is_consensus': int(cds.is_consensus),
        }
        db_cds_dict[cds.id] = cds_dict

    hash_to_cds_id_map = {}
    for cds_id, cds_dict in db_cds_dict.items():
        dict_to_hash = collections.OrderedDict({
            'name': str(cds_dict['name']),
            'namespace': str(cds_dict['namespace']),
            'is_consensus': str(cds_dict['is_consensus']),
            'exon_coords': []
        })
        # Append coords
        for cdsregion in (S.query(gbm.CdsRegion)
                          .filter(gbm.CdsRegion.cds_id == cds_id)
                          .order_by(gbm.CdsRegion.cds_order)
                          .all()):
            coord_dict = collections.OrderedDict({
                'chromosome_name': chromosome.name,
                'start': cdsregion.coordinates.start_end.lower,
                'end': cdsregion.coordinates.start_end.upper,
                'strand': cdsregion.coordinates.strand,
                'phase_start': cdsregion.phase_start,
                'phase_end': cdsregion.phase_end,
                'cds_order': cdsregion.cds_order
            })
            coord_dict = reg_dict(coord_dict)
            dict_to_hash['exon_coords'].append(coord_dict)
        # Val
        val = _hash_dict(dict_to_hash)
        if val in hash_to_cds_id_map:
            raise Exception('Unique key failed')
        hash_to_cds_id_map[val] = cds_id
    
    log.info('hash_to_cds_id_map has %d keys', len(hash_to_cds_id_map))
    return hash_to_cds_id_map

def create_hash_to_cds_accession_map(chromosome, gtf_cds_file, gtf_cdsregion_file):

    gtf_cds_dict = {}
    with open(gtf_cds_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for row_dict in reader:
            if not row_dict['chromosome_name'] == chromosome.name:
                continue
            cds_dict = {
                'name': row_dict['name'],
                'namespace': row_dict['namespace'],
                'is_consensus': int(row_dict['is_consensus'])
            }
            gtf_cds_dict[row_dict['accession']] = cds_dict 
    log.info('Found %d cds accessions', len(gtf_cds_dict))

    group = {}
    with open(gtf_cdsregion_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for row_dict in reader:
            if not row_dict['chromosome_name'] == chromosome.name:
                continue
            cdsregion_dict = collections.OrderedDict({
                'chromosome_name': row_dict['chromosome_name'],
                'start': row_dict['start'],
                'end': row_dict['end'],
                'strand': row_dict['strand'],
                'phase_start': row_dict['phase_start'],
                'phase_end': row_dict['phase_end'],
                'cds_order': int(row_dict['cds_order'])
            })
            if not row_dict['cds_accession'] in group:
                group[row_dict['cds_accession']] = []
            group[row_dict['cds_accession']].append(cdsregion_dict)
    log.info('Found %d cds groups', len(group))

    hash_to_cds_accession_map = {}
    for cds_accession, cds_dict in gtf_cds_dict.items():
        dict_to_hash = collections.OrderedDict({
            'name': str(cds_dict['name']),
            'namespace': str(cds_dict['namespace']),
            'is_consensus': str(cds_dict['is_consensus']),
            'exon_coords': []
        })
        # Append coords
        for cds_dict in sorted(group[cds_accession], key=itemgetter('cds_order')):
            coord_dict = collections.OrderedDict({
                'chromosome_name': cds_dict['chromosome_name'],
                'start': cds_dict['start'],
                'end': cds_dict['end'],
                'strand': cds_dict['strand'],
                'phase_start': cds_dict['phase_start'],
                'phase_end': cds_dict['phase_end'],
                'cds_order': cds_dict['cds_order']
            })
            coord_dict = reg_dict(coord_dict)
            dict_to_hash['exon_coords'].append(coord_dict)
        # Val
        val = _hash_dict(dict_to_hash)
        hash_to_cds_accession_map[val] = cds_accession
    
    log.info('hash_to_cds_accession_map has %d keys', len(hash_to_cds_accession_map))
    return hash_to_cds_accession_map


def compare_cdsregion(chromosome, cds_map, input_file):

    # Connect the database
    S = session(addresses.genomebrowser.database.name)

    # Fetch the cdsregions from the database for the given chromosome
    db_cdsregions_dict = {}
    for cdsregion in (S.query(gbm.CdsRegion)
                      .join(gbm.Coordinates)
                      .join(gbm.Chromosome)
                      .filter(gbm.Chromosome.id == chromosome.id)
                      .all()):
        if not cdsregion.cds_id in cds_map:
            continue
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

    compare_lists(db_cdsregions_dict, gtf_cdsregions_dict, 'cdsregion')


def compare_cds(chromosome, cds_map, input_file):
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

    # Fetch the cds from the database for given chromosome
    db_cds_dict = {}
    for cds in (S.query(gbm.CDS)
                .join(gbm.Coordinates)
                .join(gbm.Chromosome)
                .filter(gbm.Chromosome.id == chromosome.id)
                .all()):
        if not cds.id in cds_map:
            continue
        cds_dict = {
            'accession': cds_map[cds.id],
            # c. 200 of the gene ids are wrong in the db, so exclude for now
            #'gene_accession': genes_revmap[cds.gene_id],
            'name': cds.name,
            'namespace': cds.namespace,
            'is_consensus': int(cds.is_consensus),
            'chromosome_name': chromosome.name,
            'start': cds.coordinates.start_end.lower,
            'end': cds.coordinates.start_end.upper,
            'strand': cds.coordinates.strand,
        }
        db_cds_dict[cds_map[cds.id]] = cds_dict
        
    log.info('Read %d cds from database', len(db_cds_dict))

    # Read the CSV file
    gtf_cds_dict = {}
    with open(input_file, 'r') as input_fh:
        reader = csv.DictReader(input_fh)
        for cds_dict in reader:
            if cds_dict['chromosome_name'] == chromosome.name:                
                # Don't bother comparing
                cds_dict.pop('status')
                # c. 200 of the gene ids are wrong in the db, so exclude for now
                cds_dict.pop('gene_accession')
                gtf_cds_dict[cds_dict['accession']] = cds_dict
        
    log.info('Read %d cds from %s', len(gtf_cds_dict), input_file)

    compare_lists(db_cds_dict, gtf_cds_dict, 'cds')
    

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
                if int(exon_dict['phase_db_overwrite']):
                    exon_dict['phase_start'] = exon_dict['coding_start']
                    exon_dict['phase_end'] = exon_dict['coding_end']
                exon_dict.pop('coding_start')
                exon_dict.pop('coding_end')
                exon_dict.pop('phase_db_overwrite')
                tup = (exon_dict['transcript_accession'],
                       int(exon_dict['transcript_order']))
                gtf_exon_dict[tup] = exon_dict
        
    log.info('Read %d exons from %s', len(gtf_exon_dict), input_file)

    with open('gtf.out', 'w') as output_fh:
        for key in sorted(gtf_exon_dict.keys()):
            obj = gtf_exon_dict[key]
            output_fh.write(" ".join((obj['transcript_accession'],
                                      obj['transcript_order'],
                                      obj['start'],
                                      obj['end'],
                                      obj['phase_start'],
                                      obj['phase_end'])) + '\n')
    with open('db.out', 'w') as output_fh:
        for key in sorted(db_exon_dict.keys()):
            obj = db_exon_dict[key]
            output_fh.write(" ".join((obj['transcript_accession'],
                                      obj['transcript_order'],
                                      obj['start'],
                                      obj['end'],
                                      obj['phase_start'],
                                      obj['phase_end'])) + '\n')
                                       
    compare_lists(db_exon_dict, gtf_exon_dict, 'exon')




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
#        if chromosome.name != 'chr21':
#            continue
        log.info('===== Comparisons for %s =====', chromosome.name)

        # Faffing to cope with lack of unique key in cds table
        hash_to_cds_id_map = create_hash_to_cds_id_map(chromosome)
        cds_id_to_hash_map = {}
        for hash_val, cds_id in hash_to_cds_id_map.items():
            cds_id_to_hash_map[cds_id] = hash_val
        log.info('cds_id_to_hash_map has %d keys', len(cds_id_to_hash_map))
        hash_to_cds_accession_map = create_hash_to_cds_accession_map(chromosome, cds_file, cdsregions_file)
        cds_map = {}

#        pprint("===== cds_id_to_hash =====")
#        for key, obj in cds_id_to_hash_map.items():
#            pprint("%s %s" % (key, obj))
#        pprint("===== hash_to_cds_accession =====")
#        for key, obj in hash_to_cds_accession_map.items():
#            pprint("%s %s" % (key, obj))
#        exit(-1)

        for cds_id, hash_val in cds_id_to_hash_map.items():
            if hash_val in hash_to_cds_accession_map:
                cds_map[cds_id] = hash_to_cds_accession_map[hash_val]
        log.info('cds_map has %d keys', len(cds_map))

        compare_genes(chromosome, genes_file)
        compare_transcripts(chromosome, transcripts_file)
        compare_exons(chromosome, exons_file)
        compare_cds(chromosome, cds_map, cds_file)
        compare_cdsregion(chromosome, cds_map, cdsregions_file)
    
if __name__ == "__main__":
    main()
