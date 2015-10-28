"""
Test conversion of parsed snapgene object to DeskGen format.
"""

from dgparse.snapgene import parse

def test_conversion(utils):
    with utils.fixture_data("testdata/pDONR223 empty vector.dna") as f:
        parsed_ = parse(f)
    assert len(parsed_['sequence']['bases']) == 5005
