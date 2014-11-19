import sys

"""
Describe your test, its assumptions, and expected results
"""

# importing required modules (this is very ugly)
import pytest
sys.path.append('/home/rebecca/Documents/repos/dgparse/dgparse')
from snapgene import snapgene

# put some fixtures here

# Full dataset
# test datasets
@pytest.fixture(scope="module")
def fulldata():
    with open("testdata/pDONR223 empty vector.dna", "rb") as f: 
        fulldata = snapgene(f)
    return fulldata

    #return smtplib.SMTP("merlinux.eu")

with open("testdata/sequence.txt", "r") as seq:
    fullSequence = seq.read() 

# parameters of test datasets
len_DNA = 5005

# Test datasets missing specific sections to test parser error handling
#noFeatures = snapgene("testdata/test_no10.dna")
#noPrimers = snapgene("testdata/test_no5.dna")
#noOther = snapgene("testdata/test_no6.dna")
#noNotes = snapgene("testdata/test_no8.dna")


# put some tests here
class TestSnapgeneDNA:
    def testDNAMissing(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/test_no0.dna", "rb") as f:
                noDNA = snapgene(f) 
        assert excinfo.value.message == "No DNA Sequence Provided!"

    def testDNAPresent(self, fulldata):
        assert fulldata.data["DNA"] != None
        
    def testParsedDNALength(self, fulldata):
        assert len(fulldata.data["DNA"]["sequence"]) == len_DNA

    def testParsedDNASequence(self, fulldata):
        assert fulldata.data["DNA"]["sequence"] == fullSequence

    def testMethylation(self, fulldata):
        assert (fulldata.data["DNA"]["Dam"] == True) and (fulldata.data["DNA"]["Dcm"] == True) and (fulldata.data["DNA"]["EcoK1"] == True)

    def testStrandTopology(self, fulldata):
        assert (fulldata.data["DNA"]["topology"] == "circular") and (fulldata.data["DNA"]["strandedness"] == "double")

    def testTruncatedSequence(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/truncated.dna", "rb") as f:
                truncatedDNA = snapgene(f) 
        assert excinfo.value.message == "Badly formed segment or missing segment. Current segment: 29 Previous Segment: 0" 

class TestSnapgeneDescriptor:
    def testDescriptorInfo(self, fulldata):
        assert fulldata.data["descriptor"]["f_type"] == "DNA"

    def testMissingDescriptor(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/test_no9.dna", "rb") as f:
                noDescriptor = snapgene(f) 
        assert excinfo.value.message == "No snapgene Descriptor. Is this a snapgene .dna file?"

class TestOtherFeatures:
    def testFeaturesMissing(self):
        with open("testdata/test_no10.dna", "rb") as f:
            noFeatures = snapgene(f) 
        assert noFeatures.data["features"] == None

    def testPrimersMissing(self):
        with open("testdata/test_no5.dna", "rb") as f:
            noPrimers = snapgene(f) 
        assert noPrimers.data["primers"] == None

    def testOtherMissing(self):
        with open("testdata/test_no8.dna", "rb") as f:
            noOther = snapgene(f) 
        assert noOther.data["otherProperties"] == None

    def testNotesMissing(self):
        with open("testdata/test_no6.dna", "rb") as f:
            noNotes = snapgene(f) 
        assert noNotes.data["notes"] == None







#noFeatures = snapgene()
#noPrimers = snapgene("testdata/test_no5.dna")
#noOther = snapgene("testdata/test_no6.dna")
#noNotes = snapgene("testdata/test_no8.dna")
    