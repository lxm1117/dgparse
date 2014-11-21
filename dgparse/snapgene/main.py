#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

""" Parse SnapGene .dna files. See Docmentation folder for details on the file
format.
"""

# import relevant modules
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import logging
import click
import sys
import struct
import binascii
import json
from xml.dom.minidom import parseString
# from xml.dom import minidom

# setting up a dictionary of parsing functions for each segment
# if not implemented, returns the original section unparsed

from .segments import (
    decode,
    parseDNA,
    noop,
    parsePrimers,
    parseNotes,
    parseProperties,
    parseDescriptor,
    parseFeatures,
)

SEGMENT_PARSERS = {
    0: parseDNA,
    2: noop,
    3: noop,
    5: parsePrimers,
    6: parseNotes,
    8: parseProperties,
    0: parseDNA,
    9: parseDescriptor,
    10: parseFeatures,
    13: noop,
}

SEGMENT_NAME = {
    0: 'DNA',
    5: 'primers',
    8: 'otherProperties',
    6: 'notes',
    9: 'descriptor',
    10: 'features',
}


def snapgene(f):
    """
    snapgene class holds parsed Snapgene data.

    Arguments:
    f: a snapgene .dna file, opened for reading

    Attributes:
    data (dict): information parsed from the snapgene file.
    unknown (dict): unknown segments from the snapgene file (not parsed).
    Key: segment number; Value: Unparsed segment data.

    More information about the data attribute:
    data has six keys (str): DNA, descriptor, features, primers,
        otherProperties and notes, which correspond to the six file segments
        described in the SnapGene specifications. Each of these values is a
        further dictionary of key:value pairs, except for features, which is a
        list of feature dictionaries.

    Keys for the DNA dictionary:
    "sequence": DNA sequence
    "topology": DNA topology ("circular" or "linear")
    "strandedness": "single" or "double"
    "EcoK1": EcoK1 methylation state (True or False)
    "Dam": Dam methylation state (True or False)
    "Dcm": Dcm menthylation state (True or False)

    Keys for the descriptor dictionary:
    "name": Descriptor name ("Snapgene)
    "f_type": data type ("DNA" or "unknown")
    "import_version": import version (int)
    "export_version": export version (int)

    Keys for the primers dictionary:
    "primers": list of dictionaries containing primer information
        keys: description, sequence, meltingTemperature, name, simplified,
        location, annealedBases, boundStrand, recentID
    "HybridizationParams": dictionary of Hybridization Parameters
        keys: minMeltingTemperature, minContinuousMatchLen, allowMismatch)

    Keys for the notes dictionary: Description, CreatedBy, LastModified,
        ConfirmedExperimentally, references, created, CustomMapLabel, Organism,
        UseCustomMapLabel

    Keys for the otherProperties dictionary: UpstreamModification,
        UpstreamStickiness, DownstreamModification, DownstreamStickiness

    Keys for the feature dictionaries in the features list:
    "type": feature type (str)
    "name": feature name (str)
    "recentID" ID (int)
    "consecutiveTranslationNumbering": bool
    "swappedSegmentNumbering": bool
    "allowSegmentOverlaps": bool
    "readingFrame": reading frame (int)
    "Segment": sement information (dict). Keys: "color", "range", "type"



    Example usage -- DNA sequence and properties:
    mySnapgene = snapgene(myfile)
    myDNASequence = mySnapgene.data["DNA"]["sequence"]
    myDNATopology = mySnapgene.data["DNA"]["topology"]

    Example usage -- SnapGene file type and version:
    mySnapgene = snapgene(myfile)
    myFileType = mySnapgene.data["descriptor"]["f_type"]
    myFileVesrsion = mySnapgene.data["descriptor"]["import_version"]

    """
    container = {
        "DNA": None,
        "descriptor": None,
        "features": None,
        "primers": None,
        "otherProperties": None,
        "notes": None
    }
    unknown = {}

    segment = f.read(5)
    seg = None
    lastSeg = None
    while segment:
        if seg is not None:
            lastSeg = seg
        seg, seg_len = struct.unpack('>BI', segment)
        try:
            data = f.read(seg_len)
            parsedData = decode(seg, data, SEGMENT_PARSERS)
        except:
            # this is not totally robust: the error is raised because of a
            # key error following the wrong number of bytes in the previous
            # segment. It is possible that this still gives a valid key,
            # despite earlier error.
            raise Exception("Badly formed segment or missing segment. Current segment: %s Previous Segment: %s" % (seg, lastSeg))

        if seg in SEGMENT_NAME:
            if container[SEGMENT_NAME[seg]] is not None:
                errMsg = "Duplicate segments. Current segment: %s Previous segment: %s" % (seg, lastSeg)
                raise Exception(errMsg)
            else:
                container[SEGMENT_NAME[seg]] = parsedData
        else:
            # if we don't know how to parse it keep in "unknown" dictionary
            unknown[seg] = parsedData
        segment = f.read(5)

    if container["descriptor"] is None:
        raise Exception("No snapgene Descriptor. Is this a snapgene .dna file?")

    if container["DNA"] is None:
        raise Exception("No DNA Sequence Provided!")
    return container


def main():
    # TODO
    # XML parsing sections parse only those sections included in the file
    # parser(and represented in the example file), so if there are other xml
    # tags / attributes in other SnapGene files that are not represented here,
    # they will not be included.

    # Currently derived attributes from only the example file.
    # Search other .dna files for other attributes / features not included.

    # Error handling is not great:
    # custom exceptions could hide stacktrace from debugging.

    # Test output to json -- is this correct json?

    import argparse
    parser = argparse.ArgumentParser(description='Parser for SnapGene .dna \
        files. Usage: python snapgene.py /path/to/my/snapgene.dna')
    parser.add_argument('SnapGeneFile', metavar='SnapGeneFile',
                        type=str, help='filepath to SnapGene .dna file.')

    args = parser.parse_args()

    with open(args.SnapGeneFile, "r") as f:
        mySnapgene = snapgene(f)

    print json.dumps(mySnapgene, sort_keys=True, indent=4,
                     separators=(',', ': '))
