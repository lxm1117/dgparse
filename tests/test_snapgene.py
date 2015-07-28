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

# importing required modules
import os
import pytest
from dgparse.snapgene.main import parse_snapgene


@pytest.fixture(scope="module")
def parsed(utils):
    with utils.fixture_data("testdata/pDONR223 empty vector.dna") as f:
        parsed_ = parse_snapgene(f)
    return parsed_


@pytest.fixture(scope="module")
def full_sequence(utils):
    with utils.fixture_data("testdata/sequence.txt") as f:
        sequence = f.read()
    return sequence


def testDNAMissing(utils):
    with pytest.raises(Exception) as excinfo:
        with utils.fixture_data("testdata/test_no0.dna") as f:
            parse_snapgene(f)
    assert excinfo.value.message == "No DNA Sequence Provided!"


def testDNAPresent(parsed):
    assert parsed["DNA"] is not None


def testParsedDNALength(parsed):
    assert len(parsed["DNA"]["sequence"]) == 5005


def testParsedDNASequence(parsed, full_sequence):
    assert parsed["DNA"]["sequence"] == full_sequence


def testMethylation(parsed):
    assert (parsed["DNA"]["Dam"] is True) and (parsed["DNA"]["Dcm"] is True) and (parsed["DNA"]["EcoK1"] is True)


def testStrandTopology(parsed):
    assert (parsed["DNA"]["topology"] == "circular") and (parsed["DNA"]["strandedness"] == "double")


def testTruncatedSequence(utils):
    with pytest.raises(Exception) as excinfo:
        with utils.fixture_data("testdata/truncated.dna") as f:
            parse_snapgene(f)
    assert excinfo.value.message == "Badly formed segment or missing segment. Current segment: 29 Previous Segment: 0"


def testDuplicateSequence(utils):
    with pytest.raises(Exception) as excinfo:
        with utils.fixture_data("testdata/duplicate.dna") as f:
            parse_snapgene(f)
    assert excinfo.value.message == "Duplicate segments. Current segment: 0 Previous segment: 0"


"""
Test parsing of Descriptor.

f_type should be "DNA" when correct .dan file supplied. If there is no \
descriptor segemnt, raise an error
"""


def testDescriptorInfo(parsed):
    assert parsed["descriptor"]["f_type"] == "DNA"


def testMissingDescriptor(utils):
    with pytest.raises(Exception) as excinfo:
        with utils.fixture_data("testdata/test_no9.dna") as f:
            parse_snapgene(f)
    assert excinfo.value.message == "No snapgene Descriptor. Is this a snapgene .dna file?"


"""Test correct handling of missing non-essential segments"""
def testFeaturesMissing(utils):
    with utils.fixture_data("testdata/test_no10.dna") as f:
        noFeatures = parse_snapgene(f)
    assert noFeatures["features"] is None


def testPrimersMissing(utils):
    with utils.fixture_data("testdata/test_no5.dna") as f:
        noPrimers = parse_snapgene(f)
    assert noPrimers["primers"] is None


def testOtherMissing(utils):
    with utils.fixture_data("testdata/test_no8.dna") as f:
        noOther = parse_snapgene(f)
    assert noOther["other_properties"] is None


def testNotesMissing(utils):
    with utils.fixture_data("testdata/test_no6.dna") as f:
        noNotes = parse_snapgene(f)
    assert noNotes["notes"] is None
