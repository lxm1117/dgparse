# Copyright (c) 2014, Rebecca R. Murphy
# All rights reserved.

# * Redistribution and use in source and binary forms, with or without
#   modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.


import sys

"""
Describe your test, its assumptions, and expected results
"""

# importing required modules
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
class TestSnapgeneDNA(object):
    """
    Test correct parsing of DNA segment.

    When the correct .dna file is supplied, the sequence length and identity
    should match "len_DNA" and "fullSequence, the topology is circular, the
    strandedness is double and the Dam, Dcm and EcoK1 methylation states are
    all True.

    When the DNA segment is missing, an error should be raised.

    When the DNA segment is duplicated an error should be raised.

    When the segment length does not match the length described in the segment
    header, the next segment key will be unknown, so an error should be raised.

    """

    def testDNAMissing(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/test_no0.dna", "rb") as f:
                noDNA = snapgene(f)
        assert excinfo.value.message == "No DNA Sequence Provided!"

    def testDNAPresent(self, fulldata):
        assert fulldata.data["DNA"] is not None

    def testParsedDNALength(self, fulldata):
        assert len(fulldata.data["DNA"]["sequence"]) == len_DNA

    def testParsedDNASequence(self, fulldata):
        assert fulldata.data["DNA"]["sequence"] == fullSequence

    def testMethylation(self, fulldata):
        assert (fulldata.data["DNA"]["Dam"] is True) and (fulldata.data["DNA"]["Dcm"] is True) and (fulldata.data["DNA"]["EcoK1"] is True)

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


class TestSnapgeneDescriptor(object):
    """
    Test parsing of Descriptor.

    f_type should be "DNA" when correct .dan file supplied. If there is no \
    descriptor segemnt, raise an error
    """

    def testDescriptorInfo(self, fulldata):
        assert fulldata.data["descriptor"]["f_type"] == "DNA"

    def testMissingDescriptor(self):
        with pytest.raises(Exception) as excinfo:
            with open("testdata/test_no9.dna", "rb") as f:
                noDescriptor = snapgene(f)
        assert excinfo.value.message == "No snapgene Descriptor. Is this a snapgene .dna file?"


class TestOtherFeatures(object):
    """Test correct handling of missing non-essential segments"""
    def testFeaturesMissing(self):
        with open("testdata/test_no10.dna", "rb") as f:
            noFeatures = snapgene(f)
        assert noFeatures.data["features"] is None

    def testPrimersMissing(self):
        with open("testdata/test_no5.dna", "rb") as f:
            noPrimers = snapgene(f)
        assert noPrimers.data["primers"] is None

    def testOtherMissing(self):
        with open("testdata/test_no8.dna", "rb") as f:
            noOther = snapgene(f)
        assert noOther.data["otherProperties"] is None

    def testNotesMissing(self):
        with open("testdata/test_no6.dna", "rb") as f:
            noNotes = snapgene(f)
        assert noNotes.data["notes"] is None
