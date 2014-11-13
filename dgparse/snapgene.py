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
#decode_dict[0] = lambda x: decode_DNA(x)
#decode_dict[9] = lambda x: decode_descriptor(x)
#decode_dict[10] = lambda x: decode_features(x)
decode_dict[5] = lambda x: decode_primers(x)
decode_dict[8] = lambda x: decode_primers(x)
decode_dict[6] = lambda x: decode_primers(x)
#decode_dict[3] = lambda x: decode_primers(x)

def decode_DNA(data):
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

    print DNA
    print properties
    print DNA_properties
    return data

def decode_descriptor(data):
    descriptor_properties = {}

    # "snapgene" name
    descriptor_properties["name"] = data[:8]

    # file type: DNA or something else
    f_type, export_version, import_version = struct.unpack('>HHH', data[8:])
    print f_type, export_version, import_version
    if f_type == 1:
        descriptor_properties["f_type"] = "DNA"
    else:
        descriptor_properties["f_type"] = "unknown"
    descriptor_properties["export_version"] = export_version
    descriptor_properties["import_version"] = import_version
    print descriptor_properties

    return data

def decode_primers(data):
    primers = {}

    level_0 = {"recentID": str,
                "sequence": str,
                "description": str,
                "name": str,
                "minContinuousMatchLen": int,
                "minMeltingTemperature": int,
    }

    level_1 = {"meltingTemperature": int,
                "boundStrand": int,
                "simplified": bool,
                "location": str,
                "annealedBases": str,

    }

    features = ET.fromstring(data)
    for feature in features:
        print "level one", feature, feature.attrib
        for feat in feature:
            print "level two", feat, feat.attrib
            for f in feat:
                print "level three", f, f.attrib
                for ff in f:
                    print "level four!", ff, ff.attrib
                

def decode_features(data):
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
        print feature
        feature_dict = {}
        for item, func in top_level.iteritems():
            if item in feature.attrib:
                feature_dict[item] = func(feature.attrib[item])
        
        # delete after testing
        for att in feature.attrib.iterkeys():
            if att not in feature_dict:
                print "missing ", att

        feature_dict["Notes"] = []
        for feat in feature:
            segment_dict = {}
            #notes = []
            #genes = {}
            #print feat
            #print feat.tag
            #print feat.attrib
            if feat.tag == "Segment":
                for item, func in second_level.iteritems():
                    if item in feat.attrib:
                        #print item, "here!"
                        segment_dict[item] = func(feat.attrib[item])
                    #else:
                        #print item, feat.atrib
                feature_dict["Segment"] = segment_dict

            elif feat.tag == "Q":
                #if feat.attrib["name"] == "note":
                for f in feat:
                    print feat.tag, feat.attrib
                    print f.tag
                    print f.attrib
                    for k, v in f.attrib.iteritems():
                        #print k, v
                        #notes[feat.attrib["name"]] = {k:v}
                        feature_dict["Notes"].append((feat.attrib["name"], v))
            #feature_dict["Notes"] = notes               


        all_features.append(feature_dict)

    for dic in all_features:
        print dic
    return data

def decode(seg, data, parsers):
    # select correct parsing function
    func = parsers[seg]
    data = func(data)
    return data


def main():
    # TODO
    # standardize return 
    with open("../data/snapgene/pDONR223 empty vector.dna", "rb") as f:
        segment = f.read(5)
        while segment:
            seg, seg_len = struct.unpack('>BI', segment)
            #print binascii.hexlify(segment)
            print seg, seg_len
            data = f.read(seg_len)
            #print data
            snoof = decode(seg, data, decode_dict)
            #if seg == 0:
            #    print snoof
            segment = f.read(5)

if __name__ == '__main__':
    sys.exit(main())
