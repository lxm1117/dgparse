#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" Parse SnapGene .dna files. See Docmentation folder for details on the file
format.
"""
from __future__ import division
from __future__ import absolute_import
from __future__ import unicode_literals

import logging
import click
import sys
import struct
import binascii
import xml.etree.ElementTree as ET

# setting up a dictionary of parsing functions for each segment
# if not implemented, returns the original section
decode_dict = {}
for i in range(15):
    decode_dict[i] = lambda x: x
decode_dict[0] = lambda x: parseDNA(x)
decode_dict[9] = lambda x: parseDescriptor(x)
decode_dict[10] = lambda x: parseFeatures(x)
decode_dict[5] = lambda x: parsePrimers(x)
decode_dict[6] = lambda x: parseNotes(x)
decode_dict[8] = lambda x: parseProperties(x)

def parseDNA(data):
    # dictionary of DNA properties
    DNA_properties = {}

    # DNA sequence
    DNA = data[1:]
    DNA_properties["sequence"] = DNA
    
    # other properties
    properties = ord(data[0])

    # topology
    if properties & 1:
        DNA_properties["topology"] = "circular"
    else:
        DNA_properties["topology"] = "linear"

    # strandedness
    if properties & 2:
        DNA_properties["strandedness"] = "double"
    else:
        DNA_properties["strandedness"] = "single"

    # Dam methylation
    if properties & 4:
        DNA_properties["Dam"] = True
    else:
        DNA_properties["Dam"] = False

    # Dcm methylation
    if properties & 8:
        DNA_properties["Dcm"] = True
    else:
        DNA_properties["Dcm"] = False

    # EcoKI methylation
    if properties & 16:
        DNA_properties["EcoK1"] = True
    else:
        DNA_properties["EcoK1"] = False

    return DNA_properties

def parseDescriptor(data):
    descriptor_properties = {}

    # "snapgene" name
    descriptor_properties["name"] = data[:8]

    # file type: DNA or something else
    f_type, export_version, import_version = struct.unpack('>HHH', data[8:])
    #print f_type, export_version, import_version
    if f_type == 1:
        descriptor_properties["f_type"] = "DNA"
    else:
        descriptor_properties["f_type"] = "unknown"
    descriptor_properties["export_version"] = export_version
    descriptor_properties["import_version"] = import_version
    #print descriptor_properties

    return descriptor_properties

def parseNotes(data):
    level_0 = {"Synthetic": bool,
                "ConfirmedExperimentally": bool,
                "CustomMapLabel": str,
                "UseCustomMapLabel": bool,
                "Description": str,
                "Created": str,
                "CreatedBy": str,
                "LastModified": str,
                "Organism": str,
                "TransformedInot": str

    }

    ref_level = {"title": str,
                    "pubMedID": int,
                    "journal": str,
                    "authors": str

    }

    notes_dict = {}
    features = ET.fromstring(data)
    for feature in features:
        if feature.tag == "References":
            references = []
            for feat in feature:
                reference = {}
                for item, func in ref_level.iteritems():
                    if item in feat.attrib:
                        reference[item] = func(feat.attrib[item])
                references.append(reference)
            notes_dict["references"] = references
        else:
            #print feature
            for item, func in level_0.iteritems():
                if item in feature.tag:
                    notes_dict[item] = func(feature.attrib)
               
    return notes_dict

def parseProperties(data):
    level_0 = {"AdditionalSequenceProperties": str,
                "UpstreamStickiness": bool,
                    "DownstreamStickiness": bool,
                    "UpstreamModification": str,
                    "DownstreamModification": str

    }

    notes_dict = {}
    features = ET.fromstring(data)
    for feature in features:
        for item, func in level_0.iteritems():
            if item in feature.tag:
                notes_dict[item] = func(feature.attrib)
    return notes_dict

def parsePrimers(data):
    primer_dict = {}
    primers = []

    level_0 = {"recentID": str,
                "sequence": str,
                "description": str,
                "name": str,
                "minContinuousMatchLen": int,
                "minMeltingTemperature": int,
                "allowMismatch": int
    }

    level_1 = {"meltingTemperature": int,
                "boundStrand": int,
                "simplified": bool,
                "location": str,
                "annealedBases": str,
                "hybridizedRange": str,
                "BindingSite": str
    }

    features = ET.fromstring(data)
    for feature in features:
        if feature.tag == "HybridizationParams":
            # hybridization params
            feature_dict = {}
            for item, func in level_0.iteritems():
                if item in feature.attrib:
                    feature_dict[item] = func(feature.attrib[item])
            primer_dict["HybridizationParams"] = feature_dict

        if feature.tag == "Primer":
            primer = {}
            for item, func in level_0.iteritems():
                if item in feature.attrib:
                    primer[item] = func(feature.attrib[item])
                for feat in feature:
                    for item, func in level_1.iteritems():
                        if item in feat.attrib:
                            primer[item] = func(feat.attrib[item])
            primers.append(primer)
    primer_dict["primers"] = primers
    return primer_dict
                

def parseFeatures(data):
    all_features = []

    # this is hacky -- what is best structure to store parsed data?
    # dictionary would be better?
    top_level = {"name": str,
                "swappedSegmentNumbering": bool,
                "allowSegmentOverlaps": bool,
                "directionality": int,
                "consecutiveTranslationNumbering": bool,
                "readingFrame": int,
                "type": str,
                "recentID": int,
                "translationMW": float,
                "hitsStopCodon": bool
    }

    
    second_level = {"color": str, "range": str, "type": str}

    features = ET.fromstring(data)
    for feature in features:
        #print feature
        feature_dict = {}
        for item, func in top_level.iteritems():
            if item in feature.attrib:
                feature_dict[item] = func(feature.attrib[item])
        
        feature_dict["Notes"] = []
        for feat in feature:
            segment_dict = {}
            if feat.tag == "Segment":
                for item, func in second_level.iteritems():
                    if item in feat.attrib:
                        segment_dict[item] = func(feat.attrib[item])
                feature_dict["Segment"] = segment_dict

            elif feat.tag == "Q":
                for f in feat:
                    for k, v in f.attrib.iteritems():
                        feature_dict["Notes"].append((feat.attrib["name"], v))
        all_features.append(feature_dict)

    return all_features

def decode(seg, data, parsers):
    # select correct parsing function
    func = parsers[seg]
    data = func(data)
    return data

class snapgene:
    def __init__(self, fl):
        self.DNA = None
        self.descriptor = None
        self.features = None
        self.primers = None
        self.otherProperties = None
        self.notes = None
        self.unknown = []
        
        with open(fl, "rb") as f:
            segment = f.read(5)
            while segment:
                try:
                    seg, seg_len = struct.unpack('>BI', segment)
                    #print seg, seg_len
                    data = f.read(seg_len)
                    snoof = decode(seg, data, decode_dict)

                    # assign to correct attribute
                    # what is a more pythonic way to do this?
                    if seg == 0:
                        self.DNA = snoof
                    elif seg == 9:
                        self.descriptor = snoof
                    elif seg == 10:
                        self.features = snoof
                    elif seg == 5:
                        self.primers = snoof
                    elif seg == 8:
                        self.otherProperties = snoof
                    elif seg == 6:
                        self.notes = snoof
                    else:
                        self.unknown.append(snoof)
                    segment = f.read(5)
                except:
                    raise Exception("Badly formed segment or missing segment.")

            if self.descriptor == None:
                raise Exception("No snapgene Descriptor. Is this a snapgene .dna file?")

            if self.DNA == None:
                raise Exception("No DNA Sequence Provided!")



def main():
    # TODO
    # fix notes in features: make dict
    # write tests
    mySnapgene = snapgene("../data/snapgene/pDONR223 empty vector.dna")
    print mySnapgene.DNA
    print mySnapgene.descriptor
    print mySnapgene.features
    print mySnapgene.primers
    print "other:", mySnapgene.otherProperties
    print "notes:", mySnapgene.notes
    #print mySnapgene.unknown

if __name__ == '__main__':
    sys.exit(main())
