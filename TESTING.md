# Testing Protocol

## Introduction
Testing is a vital part of the software engineering process and crucial to 
having a robust parser library. This document outlines how to test parsers for
inclusion in this library.

## Methods

### The sequence hash test
Extract the sequence using a candidate parser and compute its SHA1 hash. Manually
extract that sequence or use a known-good parser on a paired file and calculate 
the hash of this "control" sequence. Assert the two hashes are equal. Assert the
lengths of the sequences are equal.

### Smoke testing
To test the robustness of a parser, a *rogues gallery* of known bad files has been
included. Each of these files has at least one problem that we've observed in 
out in the wild. Have your parser iterate over this directory and assert no 
errors are raised. 

### The matched pair import export test
The **data** directory contains matched pairs of files in common formats like
Genbank and Fasta and less common format. At a basic level a new parser should
aim to extract as much information as possible from the novel file format. Then 
if the parsed data were to be exported as a genbank file, the two files should
be reasonably similar. 