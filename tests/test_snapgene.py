from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import
import sys

"""
Describe your test, its assumptions, and expected results
"""

# importing required modules (this is very ugly)
import pytest
sys.path.append('../dgparse')
from snapgene import snapgene

# put some fixtures here

# Full dataset
fulldata = snapgene("testdata/pDONR223 empty vector.dna")
len_DNA = 5005
with open("testdata/sequence.txt", "r") as seq:
    fullSequence = seq.read() 

# Test datasets missing specific sections to test parser error handling
noDNA = "testdata/test_no0.dna"
noDescriptor = "testdata/test_no9.dna"
truncatedDNA = "testdata/truncated.dna"
noFeatures = snapgene("testdata/test_no10.dna")
noPrimers = snapgene("testdata/test_no5.dna")
noOther = snapgene("testdata/test_no6.dna")
noNotes = snapgene("testdata/test_no8.dna")


# put some tests here
class TestSnapgeneDNA:
    def testDNAMissing(self):
        with pytest.raises(Exception) as excinfo:
            snapgene(noDNA)
        assert excinfo.value.message == "No DNA Sequence Provided!"

    def testDNAPresent(self):
        assert fulldata.DNA != None
        
    def testParsedDNALength(self):
        assert len(fulldata.DNA["sequence"]) == len_DNA

    def testParsedDNASequence(self):
        assert fulldata.DNA["sequence"] == fullSequence

    def testMethylation(self):
        assert (fulldata.DNA["Dam"] == True) and (fulldata.DNA["Dcm"] == True) and (fulldata.DNA["EcoK1"] == True)

    def testStrandTopology(self):
        assert (fulldata.DNA["topology"] == "circular") and (fulldata.DNA["strandedness"] == "double")

    def testTruncatedSequence(self):
        with pytest.raises(Exception) as excinfo:
            snapgene(truncatedDNA)
        assert excinfo.value.message == "Badly formed segment or missing segment. Current segment: 29 Previous Segment: 0" 




class TestSnapgeneDescriptor:
    def testDescriptorInfo(self):
        assert fulldata.descriptor["f_type"] == "DNA"

    def testMissingDescriptor(self):
        with pytest.raises(Exception) as excinfo:
            snapgene(noDescriptor)
        assert excinfo.value.message == "No snapgene Descriptor. Is this a snapgene .dna file?"


    