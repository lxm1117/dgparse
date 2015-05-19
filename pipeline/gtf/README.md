#GTF Parser

##Introduction

##Method

##Regression test

The parser output was checked against a database provisioned with GRCh38.p2 data:

- sequences ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38.p2/Primary_Assembly/assembled_chromosomes/FASTA/
- Ensembl annotations anonymous@ensembldb.ensembl.org/homo_sapiens_core_78_38 GRCh38.p2
- CCDS data ftp://ftp.ncbi.nih.gov/pub/CCDS/current_human/CCDS.current.txt on April 18th 2015

The GTF file used was

- ftp://ftp.ensembl.org/pub/release-78/gtf/homo_sapiens/Homo_sapiens.GRCh38.78.gtf.gz

Some obvious discrepancies were:

- semi-colons in the database gene names (13 accessions) and transcripts (45 accessions) problem seems to originate at Ensembl but the GTF is not affected
- records in the database cds table linking to the wrong gene (207 cases) attributed to erroneous gene name mappings in the CCDS loader

The above differences were noted, then the comparisons were excluded from the test. 
The modified test was re-run, with results as follows:

1. identical data for the gene table
2. identical data for the transcript table
3. the parser finds a small number (~1%) of additional entries for the cds table; data is identical in the overlap set
4. the parser finds a small number (~1%) of additional entries for the cdsregion table; data is identical in the overlap set
5. the exon.end_phase for the stop codon in the exon table differs in a 0.1% of cases

Differences 3 and 4 are attributed to better capture of the relationships between coding sequence, transcript accession and CCDS identifier via the GTF file compared to the existing process involving matching names and CCDS ids from the NCBI CCDS file with Ensembl.

Difference 5 is very minor. The end_phase of the stop codon is almost always set to -1 in the database, following the convention used elsewhere for the non-coding end of an exon, but very occasionally (0.1% of exons) it is not set to -1. In contrast, the parser calculates the exon.end_phase from based on frame in the CDS file, and sets it to -1 at the end of the stop codon.

##Data

The files regression.log and regression.out in the ./tests directory provide further details and illustrate some of the differences.
