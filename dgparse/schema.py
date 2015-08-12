#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generic Schema for the DeskGen data model to support validation and
serialization.

Note: this schema is written in generic Python from a client-oriented
perspective. It is not intended as a replacement for the proprietary DTG data
model, which remains the Source of Truth.

Only attributes which are of general **END USER** interest should be shown here.
Implementation should be hidden.

Also, these Objects are essentially descriptive in nature. They don't actually
do the heavy lifting of storage. Instead a shared unmarshall (constructor) method
is called that returns a dictionary (or optionally dictionary-like) data structure.

Ultimately this module should be automatically generated from the backend
Schema.
"""
import hashlib
import json
import re
from marshmallow import Schema, fields, pre_load, validates, pre_dump

from dgparse import exc
from dgparse.sequtils import NOT_UNAMBIG_DNA, NOT_DNA, AMBIG_CHAR
# Start with the primitives and simple elements then build up


class BaseProperties(fields.Field):
    """
    Represents a key value store (JSON object) or certain properties
    """


class SequenceSchema(Schema):
    """
    An array of characters selected from a finite alphabet used
    to represent the linear structure of a biopolymer.

    Sequence is a child attribute so the underlying implementation may be
    altered provided the same interface is supported.
    """
    accession = fields.String(load_only=True)
    alphabet = fields.String(default=b'ACGT', load_only=True)  # Must be in lexographic order
    bases = fields.String()  # really a property

    @pre_load
    def get_accession(self, data):
        try:
            bases = data.get('bases').replace('\n', '').upper()
            data['accession'] = hashlib.sha1(bases).hexdigest()
            data['bases'] = bases
        except AttributeError:
            import ipdb;ipdb.set_trace()
            pass  # let the validator catch this
        finally:
            return data

    @validates('bases')
    def validate_bases(self, obj):
        if len(obj) < 12:   # no null sequences
            raise exc.NoSequence("No sequence provided.")
        hit = NOT_UNAMBIG_DNA.search(obj)
        if hit:
            msg = "Non-IUPAC Unambiguous DNA base found at {0}".format(hit.regs[0][0])
            raise exc.IllegalCharacter(msg)


class CoordinatesSchema(Schema):
    sequence_accession = fields.String()
    start = fields.Int()
    end = fields.Int  # aka lower
    strand = fields.Int


class PatternSchema(SequenceSchema):
    """
    A continuous regular pattern of bases used to define a feature
    """
    alphabet = fields.String(default=b'ACGNT', load_only=True)

    @validates('bases')
    def validate_bases(self, obj):
        if len(obj) < 4:
            raise exc.NoSequence("No Pattern Provided")
        hit = NOT_DNA.search(obj)
        if hit:
            msg = "Non-IUPAC Ambiguous DNA bases found at {0}".format(hit.regs[0][0])
            raise exc.IllegalCharacter(msg)


class RepositorySchema(Schema):
    """
    Represents a collection of biological objects either in vitro, such as
    an Inventory or in vivo such as a Genome
    """
    name = fields.String(required=True)


class BaseRepositoryItemSchema(Schema):
    """Base Inventory Item Class that defines attributes present on all
    repository items.
    The User is always assumed to be the current user writing or running the
    program.
    """
    accession = fields.String(required=True)  # A unique key for BIO objects
    created = fields.DateTime()
    modified = fields.DateTime()
    category = fields.String()  # the "Type"
    name = fields.String(required=True)
    repository = fields.Nested(RepositorySchema)  # defines accession namespace
    description = fields.String()
    properties = fields.Raw()  # General Key Value Store

    @pre_load
    def make_properties(self, data):
        """Grab out extra data and stuff it into the properties dictionary or
        it will be lost"""
        properties = data.get('properties', {})
        if isinstance(properties, basestring):
            try:
                properties = json.loads(properties)
            except ValueError:
                import ipdb;ipdb.set_trace()
                return data
        for key, val in data.iteritems():
            if val is None:
                continue # skip blank fields
            if val is "":
                continue
            if key in self.fields:
                continue  # we'll grab it properly later
            if key not in properties:  #don't over-write
                properties[key] = val
        data['properties'] = properties
        return data

    @pre_dump
    def move_to_properties(self, obj):
        """Move unsupported attributes to the properties dictionary"""
        if isinstance(obj, list):
            return
        for key in self.fields.iterkeys():
            field_obj = self.fields[key]
            if field_obj.load_only and key in obj:  # convention for things we don't serialize
                obj['properties'][key] = obj.pop(key)


# The Physical Things
class BaseRepositoryFileSchema(BaseRepositoryItemSchema):
    """
    The Source files uploaded for various items.
    """
    format_ = fields.String()  # file format NOT content type/category
    parsed = fields.Bool()
    record_type = fields.String()
    contents = fields.Raw()
    file_path = fields.String()
    size = fields.Integer()


class SequencingReadFile(BaseRepositoryFileSchema):
    """
    A file representing the raw, experimental support for the existence of a
    biopolymer.
    """
    molecule_accession = fields.String(required=True)  # pseudo fkey back to molecule


class BaseAnnotationSchema(BaseRepositoryItemSchema):
    quality = fields.Float()  # How true is this annotation
    # relationships
    order = fields.Int(default=0)
    length = fields.Int()
    alignment = fields.String()
    regions = fields.Nested('self', many=True, order=True)
    # relationships
    coordinates = fields.Nested(CoordinatesSchema)
    # Properties
    sequence = fields.Nested(SequenceSchema)


class BaseMoleculeSchema(BaseRepositoryItemSchema):
    """
    Represents the record of some sort of biologically important Molecule.
    """
    # attributes
    strand_count = fields.Int(default=1, load_only=True)
    sequence_file_path = fields.String(load_only=True)  # if the sequence is not present
    location = fields.String()
    # properties
    length = fields.Integer(required=True)
    mol_weight = fields.Float()
    concentration = fields.Float(load_from='conc_um')
    concentration_units = fields.String(default='ng/ul', load_only=True)

    # relationships
    sequence = fields.Nested(SequenceSchema, required=True)
    annotations = fields.Nested(BaseAnnotationSchema, many=True)
    reads = fields.Nested(SequencingReadFile, many=True)

    @pre_load
    def get_length(self, data):
        if 'length' in data and data['length'] > 0:
            return data
        try:
            sequence = data['sequence']
            data['length'] = len(sequence['bases'])
        except TypeError:
            data['length'] = None
        except KeyError:
            # raise exc.NoSequence("No Sequence Provided")
            data['length'] = None
        return data


class DnaMoleculeFileSchema(Schema):

    name = fields.String(default="DirectUpload")
    contents = fields.String(default='EMPTY')
    format_ = fields.String(default='fasta')
    parsed = fields.Boolean(default=False)


class DnaMoleculeSchema(BaseMoleculeSchema):
    date_stored = fields.DateTime()  # When the banking took place
    quality = fields.Float(load_only=True)
    sequencing_notes = fields.String(load_only=True)
    strand_count = fields.Integer(default=2, load_only=True)
    notebook_page = fields.String(default="Undefined", load_only=True)


class DnaPlasmidSchema(DnaMoleculeSchema):
    """
    A banked plasmid in an Inventory.
    """
    is_available = fields.Bool(load_only=True)
    is_circular = fields.Boolean(True)
    integration_locus = fields.String()
    carrier_strain = fields.String()
    dnamoleculefile = fields.Nested(DnaMoleculeFileSchema, dump_only=True)
    source_file = fields.Nested(BaseRepositoryFileSchema)
    files = fields.Nested(BaseRepositoryFileSchema, many=True)

    @pre_dump
    def make_fake_file(self, data):
        if 'dnamoleculefile' or 'dnamoleculefile_id' not in data:
            fake_file = {
                'name': "DirectUpload",
                'contents': 'EMPTY',
                'parsed': True
            }
            data['dnamoleculefile'] = fake_file
        return data


class DnaConstruct(DnaMoleculeSchema):
    """
    A linear fragment of DNA in an inventory.
    """
    is_circular = fields.Boolean(default=False)


class SequenceModSchema(Schema):
    """
    Represents some sort of modification to a sequence
    """
    position = fields.Integer()
    symbol = fields.String()


class DnaOligoSchema(DnaMoleculeSchema):
    """
    A linear, single-stranded piece of DNA
    """
    #accession must start with o
    external_identifier = fields.String(load_only=True)
    is_circular = fields.Boolean(False, load_only=True)  # for now
    strand_count = fields.Integer(default=1, load_only=True)
    concentration_units = fields.String(default=u'uM', load_only=True)
    t_melt = fields.Float(load_only=True)
    t_melt_method = fields.String(default='Nearest Neighbor', load_only=True)
    works = fields.Boolean(load_only=True)
    notebook_xref = fields.String(load_only=True)
    sequence = fields.Nested(SequenceSchema, required=True)
    modifications = fields.Nested(SequenceModSchema, many=True, load_only=True)
    delta_g = fields.Float(allow_none=True, load_only=True)
    target = fields.String(load_only=True)

    @pre_load
    def extract_modifications(self, data):
        clean_bases = []
        modifications = []
        if 'sequence' not in data:
            return data
        if not data['sequence']['bases']:
           return data  # let the validator handle it
        for i, base in enumerate(data['sequence']['bases']):
            if base not in AMBIG_CHAR:
                modifications.append({'position': i, 'symbol': base})
            else:
                clean_bases.append(base)
        data['sequence']['bases'] = ''.join(clean_bases)
        data['modifications'] = modifications
        return data

    @pre_load
    def parse_concentration(self, data):
        conc_string = data.pop('concentration', None)
        if conc_string is None:
            conc_string = data.pop('conc_um', None)
        if isinstance(conc_string, basestring):
            pattern = r'(\d+)\s?(\w+)'
            match = re.search(pattern, conc_string)
            if match:
                data['concentration'] = match.group(1)
                data['conc_um'] = match.group(1)
                data['concentration_units'] = match.group(2)
        return data

    @pre_load
    def clean_t_melt(self, data):
        if 't_melt' in data:
            t_melt = data.pop('t_melt')
            if isinstance(t_melt, basestring):
                value, units = t_melt.split(u'°')
                data['t_melt'] = value
            else:
                data['t_melt'] = t_melt
        return data

    @pre_dump
    def put_length(self, data):
        if data['length'] < 1:
            data['length'] = len(data['sequence']['bases'])
        if 'concentration' in data:
            data.pop('concentration')  # not supported
        if 'concentration_units' in data:
            data.pop('concentration_units')
        return data


class DnaPrimerSchema(DnaOligoSchema):
    """
    A DNA Oligo that is used to prime a PCR reaction.
    """
    # accession must start with "m"
    homology_sequence = fields.String(load_only=True)
    priming_sequence = fields.String(load_only=True)


class Chromosome(DnaMoleculeSchema):
    """
    Represents a Chromosome
    """
    genome = fields.String()


class GeneSchema(BaseAnnotationSchema):
    """
    Represents a gene and it's child objects.
    """
    biotype = fields.String(default="Unknown")


class TranscriptSchema(BaseAnnotationSchema):
    """
    Represents a Transcript
    """


class ExonSchema(BaseAnnotationSchema):
    """
    Represents an Exon.
    """


class NucleaseCutSiteSchema(BaseAnnotationSchema):
    """
    The target of a
    """


class PolypeptideSchema(BaseRepositoryItemSchema):
    """
    Represents a protein
    """


class NucleaseSchema(PolypeptideSchema):
    """
    A protein that cuts DNA
    """


class GuideRnaCutSchema(NucleaseCutSiteSchema):
    """
    Guide Cut Site Schema
    """
    prefix = fields.String()
    protospacer = fields.String()
    pam = fields.String()
    # sequence is the full sequence
    suffix = fields.String()
    coordinates = fields.Nested(CoordinatesSchema)
    cut_site = fields.Integer()
    centerpoint = fields.Integer()
    activity_score = fields.Float()  # The on-target activity score
    specificity_score = fields.Float()  #The aggregate off target score
    offtargets = fields.Nested('self', many=True)  # simpler
    failed_filter_name = None
    coding_mismatch = None
    noncoding_mismatch = None
    gc_region = None
    nuclease = fields.Nested(NucleaseSchema)  # change to full nuclease
    sequence = fields.Nested(SequenceSchema)


class BaseFeatureSchema(BaseRepositoryItemSchema):
    pattern = fields.Nested(PatternSchema, required=True)
    length = fields.Integer(required=True)


class DnaFeatureCategorySchema(Schema):
    name = fields.String()


class DnaFeatureSchema(BaseFeatureSchema):
    accession = fields.String(required=True)
    dnafeaturecategory_id = fields.Integer()
    pattern = fields.Nested(PatternSchema, required=True)
    category = fields.Nested(DnaFeatureCategorySchema)

    @pre_load
    def make_accession(self, obj):
        if 'accession' not in obj:
            accession = '/'.join([obj['category']['name'], obj['name']])
            obj['accession'] = accession
        return obj

    @pre_dump
    @pre_load
    def get_length(self, data):
        if 'pattern' in data and 'bases' in data['pattern']:
            pattern = data.get('pattern')
            data['length'] = len(pattern['bases'])
        return data


class BaseDesignSchema(BaseRepositoryItemSchema):
    sequence = fields.Nested(SequenceSchema)
    annotations = fields.Nested(BaseAnnotationSchema, many=True)


class DnaDesignSchema(BaseDesignSchema):
    sequence = fields.String(required=True)

    @pre_load
    def adapt_sequence(self, data):
        # TODO remove this after updating AC schema
        sequence = data.pop('sequence')
        if isinstance(sequence, dict):
            try:
                data['sequence'] = sequence.get('bases')
            except KeyError:
                data['sequence'] = sequence
        else:
            data['sequence'] = sequence
        return data


class RnaGuidedNucleaseSchema(NucleaseSchema):
    """
    The actual RNA Guided Nuclease itself.
    """


class SegmentSchema(Schema):
    """
    A segment of a cloning solution
    """

