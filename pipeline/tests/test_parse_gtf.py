from nose.tools import *
from pkg_resources import resource_filename
from gtf.parse_gtf import main

import os
import tempfile
import csv
import shutil
import filecmp


OUTER_TODAY = None


def setup():
    'Setup function that runs before every test'
    print "SETUP!"


def teardown():
    'Teardown function that runs after every test'
    print "TEAR DOWN!"


def test_main():
    'Test prepare.main'
    # Work area
    work_dir = tempfile.mkdtemp()
    # Input files
    gtf_file = \
        resource_filename(__name__, 'data/gtf/chr21.78.gz')
    genome_config_file = \
        resource_filename(__name__, 'data/gtf/genome_config.csv')
    template_file = \
        resource_filename(__name__, 'data/gtf/template.json')
    compare_dir = \
        resource_filename(__name__, 'data/gtf')
    # Run main
    main({
        '--config': genome_config_file,
        '--input': gtf_file,
        '--namespace': 'ensembl',
        '--require-ccds': True,
        '--template': template_file,
        '--transform': True,
        '--output': work_dir,
        '--log': 'INFO',
    })
    # Check files
    filenames = [ "db.{0}.csv".format(key) for key in [
        'gene', 'cds', 'cdsregion', 'transcript', 'exon']]
    for filename in filenames:
        assert filecmp.cmp(
            work_dir + '/' + filename, 
            compare_dir + '/' + filename)
    # Clean up                                                                                                                                                                                \
    shutil.rmtree(work_dir)
