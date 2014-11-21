import struct
import xml.etree.ElementTree as ET

noop = lambda x: x


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
    # repared_data = parseString(data)
    # print repared_data.toprettyxml(indent="\t")
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
            for item, func in level_0.iteritems():
                if item in feature.tag:
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
        feature_dict = {}
        for item, func in top_level.iteritems():
            if item in feature.attrib:
                feature_dict[item] = func(feature.attrib[item])

        feature_dict["Notes"] = {}
        for feat in feature:
            segment_dict = {}
            if feat.tag == "Segment":
                for item, func in second_level.iteritems():
                    if item in feat.attrib:
                        segment_dict[item] = func(feat.attrib[item])
                feature_dict["Segment"] = segment_dict

            elif (feat.tag == "Q") or (feat.tag == "Qualifier"):
                for f in feat:
                    for k, v in f.attrib.iteritems():
                        feature_dict["Notes"][feat.attrib["name"]] = v
        all_features.append(feature_dict)

    return all_features


def decode(seg, data, parsers):
    # select correct parsing function
    func = parsers[seg]
    data = func(data)
    return data
