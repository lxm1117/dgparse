import sys

"""
Describe your test, its assumptions, and expected results
"""

# importing required modules (this is very ugly)
import pytest
sys.path.append('../dgparse')
from snapgene import snapgene

# put some fixtures here

# Full dataset for evaluating 
@pytest.fixture(scope="module")
def fulldata():
    with open("testdata/pDONR223 empty vector.dna", "rb") as f: 
        fulldata = snapgene(f)
    return fulldata

# parameters of test datasets
len_DNA = 5005

with open("testdata/sequence.txt", "r") as seq:
    fullSequence = seq.read() 


# put some tests here
class TestSnapgeneDNA:
    """
    Test correct parsing of DNA segment.

    When the correct .dna file is supplied, the sequence length and identity should match "len_DNA" and "fullSequence, the topology is circular, the strandedness is double and the Dam, Dcm and EcoK1 methylation states are all True.

    When the DNA segment is missing, an error should be raised.

    When the DNA segment is duplicated an error should be raised.

    When the segment length does not match the length described in the segment header, the next segment key will be unknown, so an error should be raised.
    
    """

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

    def testDuplicateSequence(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/duplicate.dna", "rb") as f:
                truncatedDNA = snapgene(f) 
        assert excinfo.value.message == "Duplicate segments. Current segment: 0 Previous segment: 0" 
 

class TestSnapgeneDescriptor:
    """
    Test parsing of Descriptor. 

    f_type should be "DNA" when correct .dan file supplied. If there is no descriptor segemnt, raise an error
    """
    def testDescriptorInfo(self, fulldata):
        assert fulldata.data["descriptor"]["f_type"] == "DNA"

    def testMissingDescriptor(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/test_no9.dna", "rb") as f:
                noDescriptor = snapgene(f) 
        assert excinfo.value.message == "No snapgene Descriptor. Is this a snapgene .dna file?"

class TestOtherFeatures:
    """Test correct handling of missing non-essential segments"""
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
    