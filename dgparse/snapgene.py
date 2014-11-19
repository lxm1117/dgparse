#!/usr/bin/env python
# -*- coding: utf-8 -*-
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
import xml.etree.ElementTree as ET
import json
from xml.dom.minidom import parseString
#from xml.dom import minidom

# setting up a dictionary of parsing functions for each segment
# if not implemented, returns the original section unparsed

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
                "CustomMapLabel": unicode,
                "UseCustomMapLabel": bool,
                "Description": unicode,
                "Created": unicode,
                "CreatedBy": unicode,
                "LastModified": unicode,
                "Organism": unicode,
                "TransformedInot": unicode

    }

    ref_level = {"title": str,
                    "pubMedID": int,
                    "journal": str,
                    "authors": str

    }

    notes_dict = {}
    features = ET.fromstring(data)
    repared_data = parseString(data)
    #print repared_data.toprettyxml(indent="\t")
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
                    #print feature.text
                    notes_dict[item] = func(feature.text)
               
    return notes_dict

def parseProperties(data):
    level_0 = {"AdditionalSequenceProperties": unicode,
                "UpstreamStickiness": bool,
                    "DownstreamStickiness": bool,
                    "UpstreamModification": unicode,
                    "DownstreamModification": unicode

    }

    properties_dict = {}
    features = ET.fromstring(data)
    for feature in features:
        for item, func in level_0.iteritems():
            if item in feature.tag:
                properties_dict[item] = func(feature.text)
    return properties_dict

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
                        feature_dict["Notes"].append({feat.attrib["name"]: v})
        all_features.append(feature_dict)

    return all_features

def decode(seg, data, parsers):
    # select correct parsing function
    func = parsers[seg]
    data = func(data)
    return data

class snapgene:
    """
    snapgene class holds parsed Snapgene data.

    Arguments: 
    f: a snapgene .dna file, opened for reading

    Attributes:
    data (dict): information parsed from the snapgene file.
    unknown (dict): unknown segments from the snapgene file (not parsed). Key: segment number; Value: Unparsed segment data.

    More information about the data attribute:
    data has six keys (str): DNA, descriptor, features, primers, otherProperties and notes, which correspond to the six file segments described in the SnapGene specifications. Each of these values is a further dictionary of key:value pairs, except for features, which is a list of feature dictionaries.

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
    "primers": list of dictionaries containing primer information (keys: description, sequence, meltingTemperature, name, simplified, location, annealedBases, boundStrand, recentID)
    "HybridizationParams": dictionary of Hybridization Parameters (keys: minMeltingTemperature, minContinuousMatchLen, allowMismatch)

    Keys for the notes dictionary: Description, CreatedBy, LastModified, ConfirmedExperimentally, references, created, CustomMapLabel, Organism, UseCustomMapLabel

    Keys for the otherProperties dictionary: UpstreamModification, UpstreamStickiness, DownstreamModification, DownstreamStickiness

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
    map_dict = {0: "DNA",
            9: "descriptor",
            10: "features",
            5: "primers",
            8: "otherProperties",
            6: "notes"
    
    }

    decode_dict = {}
    for i in [2, 3, 13]:
        decode_dict[i] = lambda x: x
    decode_dict[0] =  parseDNA
    decode_dict[9] = parseDescriptor
    decode_dict[10] = parseFeatures
    decode_dict[5] = parsePrimers
    decode_dict[6] = parseNotes
    decode_dict[8] = parseProperties

    def __init__(self, f):
        self.data = {"DNA": None,
                        "descriptor": None,
                        "features": None,
                        "primers": None,
                        "otherProperties": None,
                        "notes": None
                        }
        self.unknown = {}
            
        segment = f.read(5)
        seg = None
        lastSeg = None
        while segment:
            if seg != None:
                lastSeg = seg
            seg, seg_len = struct.unpack('>BI', segment)
            #print seg
            try:
                data = f.read(seg_len)
                parsedData = decode(seg, data, snapgene.decode_dict)
            except:
                # this is not totally robust: the error is raised because of a key error following the wrong number of bytes in the previous segment. It is possible that this still gives a valid key, despite earlier error.
                raise Exception("Badly formed segment or missing segment. Current segment: %s Previous Segment: %s" %(seg, lastSeg))
                #raise

            if seg in snapgene.map_dict:
                if self.data[snapgene.map_dict[seg]] != None:
                    raise Exception("Duplicate segments. Current segment: %s Previous segment: %s" %(seg, lastSeg))
                else:
                    self.data[snapgene.map_dict[seg]] = parsedData
            else:
                self.unknown[seg] = parsedData
            segment = f.read(5)

        if self.data["descriptor"] == None:
            raise Exception("No snapgene Descriptor. Is this a snapgene .dna file?")

        if self.data["DNA"] == None:
            raise Exception("No DNA Sequence Provided!")

def main():
    # TODO
    # XML parsing sections parse only those sections included in the file parser (and represented in the example file), so if there are other xml tags / attributes that are not represented here, they will not be included.

    # Currently derived attributes from only the example file. Search other .dna files for other attributes / features not included.

    # Error handling is poor. Custom exception hides stacktrace from debugging.

    # Test output to json -- is this correct json?
    import argparse
    parser = argparse.ArgumentParser(description='Parser for SnapGene .dna files')
    parser.add_argument('SnapGeneFile', metavar='SnapGeneFile', type=str, help='filepath to SnapGene .dna file.')

    args = parser.parse_args()

    with open(args.SnapGeneFile, "r") as f:
        mySnapgene = snapgene(f)
    
    #print "DNA:", mySnapgene.data["DNA"]
    #print mySnapgene.data["DNA"]["sequence"]
    
    #for it in mySnapgene.data["features"]:
    #    for k, v in it.iteritems():
    #        print k, v

    print json.dumps(mySnapgene.data, sort_keys=True, indent=4, separators=(',', ': '))

if __name__ == '__main__':
    sys.exit(main())

