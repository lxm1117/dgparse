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
import json
import re
from marshmallow import Schema, fields, pre_load, validates, pre_dump

from dgparse import exc
from dgparse.sequtils import NOT_UNAMBIG_DNA, NOT_DNA, AMBIG_CHAR, compute_sha1
# Start with the primitives and simple elements then build up


class SequenceSchema(Schema):
    """
    An array of characters selected from a finite alphabet used
    to represent the linear structure of biopolymer (DNA, RNA, Protein).

    Sequence is a child attribute so the underlying implementation may be
    altered provided the same interface is supported.
    """
    sha1 = fields.String(load_only=True)
    alphabet = fields.String(default=b'ACGT', load_only=True)  # Must be in lexographic order
    bases = fields.String()  # really a property

    @pre_load
    def compute_sha1(self, data):
        """If the sha1 hash is not provided, compute it."""
        return compute_sha1(data)

    @validates('bases')
    def validate_bases(self, obj):
        """Validate dirty bases are not passed"""
        if len(obj) < 12:   # no null sequences
            raise exc.NoSequence("No sequence provided.")
        hit = NOT_UNAMBIG_DNA.search(obj)
        if hit:
            msg = "Non-IUPAC Unambiguous DNA base found at {0}".format(hit.regs[0][0])
            raise exc.IllegalCharacter(msg)


class PatternSchema(SequenceSchema):
    """
    A continuous regular pattern of bases used to define a feature of a molecule.
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


class BaseOrganisationSchema(Schema):
    """Basic Organisation Schema"""
    name = fields.String()
    domain = fields.String()
    admin_email = fields.Email()


class BaseUserSchema(Schema):
    """Basic User Schema"""
    name = fields.String()
    email = fields.Email()


class RepositorySchema(Schema):
    """
    Represents a collection of biological objects either in vitro, such as
    an Inventory or in vivo such as a Genome
    """
    name = fields.String(required=True)


class CoordinatesSchema(Schema):
    """
    Absolute Coordinates in a molecule
    """
    start_end = fields.Raw(required=True)
    strand = fields.Integer(required=True)


class BaseRepositoryItemSchema(Schema):
    """Base Inventory Item Class that defines attributes present on all
    repository items.
    The User is always assumed to be the current user writing or running the
    program.
    """
    accession = fields.String(required=True)  # A unique key for BIO objects
    created = fields.String()  # Crudpile constructor current handles these
    modified = fields.String()
    category = fields.String()  # the "Type"
    name = fields.String(required=True)
    repository = fields.Nested(RepositorySchema)  # defines sha1 namespace
    description = fields.String(required=False)
    notes = fields.String(required=False, allow_none=True)
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


class BaseRepositoryFileSchema(BaseRepositoryItemSchema):
    """
    The Source files uploaded for various items.
    """
    format_ = fields.String()  # file format NOT content type/category
    parsed = fields.Bool()
    record_type = fields.String()
    contents = fields.Raw()
    filepath = fields.String()
    size = fields.Integer()


class SequencingReadFile(BaseRepositoryFileSchema):
    """
    A file representing the raw, experimental support for the existence of a
    biopolymer.
    """
    molecule_accession = fields.String(required=True)  # pseudo fkey back to molecule


class BaseAnnotationSchema(BaseRepositoryItemSchema):
    """
    Represents an annotation of a region of a molecule, something that has
    coordinates in a molecule.
    """
    quality = fields.Float(default=0.5, allow_none=True)  # How true is this annotation
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
    location = fields.String()
    is_available = fields.Bool(default=True)
    is_circular = fields.Bool(default=False)  # *topology*
    # properties
    length = fields.Integer(required=True)
    mol_weight = fields.Float()
    concentration = fields.Float(load_from='conc_um')
    concentration_units = fields.String(default='ng/ul', load_only=True)
    volume = fields.Float()
    volume_units = fields.String(default='ul')
    # relationships
    sequence = fields.Nested(SequenceSchema, required=True)
    annotations = fields.Nested(BaseAnnotationSchema, many=True)
    reads = fields.Nested(SequencingReadFile, many=True)

    @pre_load
    def get_length(self, data):
        """
        Compute the length of the molecule
        :data:
        """
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
    """An accessory data file for a DNA Molecule"""
    name = fields.String(default="DirectUpload")
    contents = fields.String(default='EMPTY')
    format_ = fields.String(default='fasta')
    parsed = fields.Boolean(default=False)


class DnaMoleculeSchema(BaseMoleculeSchema):
    """
    A double stranded DNA molecule.
    """
    date_stored = fields.DateTime(allow_none=True)  # When the banking took place
    quality = fields.Float(load_only=True, default=0.5)
    sequencing_notes = fields.String(load_only=True)
    strand_count = fields.Integer(default=2, load_only=True)
    notebook_page = fields.String(default="Undefined", load_only=True, allow_none=True)
    location = fields.String(allow_none=True)
    description = fields.String(allow_none=True)


class DnaPlasmidSchema(DnaMoleculeSchema):
    """
    A banked plasmid in an Inventory.
    """
    is_available = fields.Bool(load_only=True)
    category = fields.String(default='plasmid')
    is_circular = fields.Boolean(True)
    integration_locus = fields.String()
    carrier_strain = fields.String()
    dnamoleculefile = fields.Nested(DnaMoleculeFileSchema, dump_only=True)
    source_file = fields.Nested(BaseRepositoryFileSchema)
    files = fields.Nested(BaseRepositoryFileSchema, many=True)

    @pre_dump
    def make_fake_file(self, data):
        """
        A temp work around as files are currently required of all molecules
        in the database.
        """
        if 'dnamoleculefile' or 'dnamoleculefile_id' not in data:
            fake_file = {
                'name': "DirectUpload",
                'contents': 'EMPTY',
                'parsed': True
            }
            data['dnamoleculefile'] = fake_file
        return data


class DnaConstructSchema(DnaMoleculeSchema):
    """
    A linear fragment of DNA in an inventory.
    """
    is_circular = fields.Boolean(default=False)
    category = fields.String(default='construct')


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
    #sha1 must start with o
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
        """
        Extract chemical modifications to the bonds or bases of a sequence.
        """
        clean_bases = []
        modifications = []
        if 'sequence' not in data:
            return data
        if not isinstance(data['sequence'], dict):
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
        """
        Figure out which sort of concentration units have been provided
        """
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
        """Parse and clean up the Melting Temperature and handle unicode"""
        if 't_melt' in data:
            t_melt = data.pop('t_melt')
            if isinstance(t_melt, basestring):
                value = t_melt.split(u'°')[0]
                data['t_melt'] = float(value)
            else:
                data['t_melt'] = t_melt
        return data

    @pre_dump
    def put_length(self, data):
        """Add the length to the oligo"""
        if 'length' not in data or data['length'] < 1:
            if 'sequence' in data:
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
    is_circular = fields.Boolean(False)
    strand_count = fields.Integer(default=1)
    concentration_units = fields.String(default=u'uM')
    t_melt_method = fields.String(default='Primer3')
    works = fields.Boolean()
    notebook_xref = fields.String()
    sequence = fields.Nested(SequenceSchema)
    modifications = fields.Nested(SequenceModSchema, many=True)
    delta_g = fields.Float()
    target = fields.String()
    # sha1 must start with "m"
    homology_sequence = fields.String(load_only=True)
    priming_sequence = fields.String(load_only=True)
    category = fields.String(default='primer')


class Chromosome(DnaMoleculeSchema):
    """
    Represents a Chromosome
    """
    category = fields.String(default='chromosome')
    genome = fields.String()


class GeneSchema(BaseAnnotationSchema):
    """
    Represents a gene and it's child objects.
    """
    biotype = fields.String(default="Unknown")  # category?


class TranscriptSchema(BaseAnnotationSchema):
    """
    Represents a Transcript
    """


class ExonSchema(BaseAnnotationSchema):
    """
    Represents an Exon.
    """


class PolypeptideSchema(BaseRepositoryItemSchema):
    """
    Represents a protein
    """


class NucleaseSchema(PolypeptideSchema):
    """
    A protein that cuts DNA
    """


class BaseCutSiteSchema(BaseRepositoryItemSchema):
    """
    The target of a nuclease is the cut site
    """
    # name, sha1, etc all included
    molecule_accession = fields.String()
    coordinates = fields.Nested(CoordinatesSchema)
    topcut = fields.Integer()  # absolute coordinates
    btmcut = fields.Integer()  # absolute coordinates
    nuclease = fields.Nested(NucleaseSchema)


class GuideRnaCutSchema(BaseCutSiteSchema):
    """
    Guide Cut Site Schema
    """
    prefix = fields.String()
    protospacer = fields.String()
    pam = fields.String()
    # sequence is the full sequence
    suffix = fields.String()
    cut_site = fields.Integer()  # absolute coordinates
    centerpoint = fields.Integer()  # absolute coordinates
    activity_score = fields.Float()  # The on-target activity score
    specificity_score = fields.Float()  #The aggregate off target score
    offtargets = fields.Nested('self', many=True)  # simpler
    failed_filter_name = None
    coding_mismatch = None
    noncoding_mismatch = None
    gc_region = None


class BaseFeatureSchema(BaseRepositoryItemSchema):
    """
    The basic definition of a biological part, a functional unit of a molecule.
    """
    pattern = fields.Nested(PatternSchema, required=True)
    length = fields.Integer(required=True)


class DnaFeatureCategorySchema(Schema):
    """
    Represents a 'soft class' of DNA Feature with a particular function.
    """
    name = fields.String()


class DnaFeatureSchema(BaseFeatureSchema):
    """
    A Dna Part Definition. A grammatical unit of DNA which has some biological
    significance.
    """
    accession = fields.String(required=True)
    dnafeaturecategory_id = fields.Integer()
    pattern = fields.Nested(PatternSchema, required=True)
    category = fields.Nested(DnaFeatureCategorySchema)

    @pre_load
    def make_accession(self, obj):
        """
        Construct a unique identifier for the feature if not provided.
        """
        if 'sha1' not in obj:
            accession = '/'.join([obj['category']['name'], obj['name']])
            obj['sha1'] = accession
        return obj

    @pre_dump
    @pre_load
    def get_length(self, data):
        """
        Construct the length of the feature from it's pattern if not provided.
        """
        if 'pattern' in data and 'bases' in data['pattern']:
            pattern = data.get('pattern')
            data['length'] = len(pattern['bases'])
        return data


class BaseDesignSchema(BaseRepositoryItemSchema):
    """
    A category of DNAMolecule that has not yet by physically constructed.
    These are used to represent the goal molecule of a genome editing or cloning
    experiment.
    """
    sequence = fields.Nested(SequenceSchema)
    annotations = fields.Nested(BaseAnnotationSchema, many=True)


class DnaDesignSchema(BaseDesignSchema):
    """
    A goal DNA molecule; something the user is trying to obtain.
    """
    sequence = fields.String(required=True)

    @pre_load
    def adapt_sequence(self, data):
        """
        Temporary adapter to put sequence on the primary object rather than
        a child table.
        """
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


class BaseRestrictionEnzymeSchema(NucleaseSchema):
    """
    A Type II restriction endonuclease enzyme
    """


class BaseExperimentItemSchema(BaseRepositoryItemSchema):
    """
    A step in a magical laboratory quest.
    """
    objects = fields.Raw()


class BaseExperimentSchema(BaseRepositoryItemSchema):
    """
    A magical laboratory quest.
    """
    items = fields.Nested(BaseExperimentItemSchema, many=True)
