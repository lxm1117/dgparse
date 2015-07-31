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
from marshmallow import Schema, fields

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
    accession = fields.String()  # The Sha1 keep a consistent interface
    alphabet = fields.String(default=b'ACGT')  # Must be in lexographic order
    bases = fields.String()  # really a property


class CoordinatesSchema(Schema):
    sequence_accession = fields.String()
    start = fields.Int()
    end = fields.Int  # aka lower
    strand = fields.Int


class PatternSchema(SequenceSchema):
    """
    A continuous regular pattern of bases used to define a feature
    """
    alphabet = fields.String(default=b'ACGNT')


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
    notes = fields.String()
    properties = fields.Raw()  # General Key Value Store


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
    molecule_accession = fields.String()  # pseudo fkey back to molecule


class BaseMoleculeSchema(BaseRepositoryItemSchema):
    """
    Represents the record of some sort of biologically important Molecule.
    """
    # attributes
    is_available = fields.Bool(default=True)
    is_circular = fields.Bool(default=False)  # *topology*
    strand_count = fields.Int(default=1)
    sequence_file_path = fields.String()  # if the sequence is not present
    location = fields.String()
    # properties
    length = fields.Int()
    mol_weight = fields.Float()
    concentration = fields.Float()
    concentration_units = fields.String(default='ng/ul')

    # relationships
    sequence = fields.Nested(SequenceSchema)
    annotations = fields.Nested(many=True)
    reads = fields.Nested(SequencingReadFile, many=True)
    source_file = fields.Nested(BaseRepositoryFileSchema)
    files = fields.Nested(BaseRepositoryFileSchema, many=True)


class DnaMoleculeSchema(BaseMoleculeSchema):

    date_stored = fields.DateTime()  # When the banking took place
    quality = fields.Float()
    sequencing_notes = fields.String()
    strand_count = fields.Integer(default=2)
    notebook_page = fields.String(default="Undefined")


class DnaPlasmidSchema(DnaMoleculeSchema):
    """
    A banked plasmid in an Inventory.
    """
    is_circular = fields.Boolean(True)
    integration_locus = fields.String()
    carrier_strain = fields.String()


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


class DnaPrimer(DnaOligoSchema):
    """
    A DNA Oligo that is used to prime a PCR reaction.
    """
    # accession must start with "m"
    homology_sequence = fields.String()
    priming_sequence = fields.String()


class Chromosome(DnaMoleculeSchema):
    """
    Represents a Chromosome
    """
    genome = fields.String()


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
    pattern = fields.Nested(PatternSchema)


class BaseDesignSchema(BaseRepositoryItemSchema):
    sequence = fields.Nested(SequenceSchema)
    annotations = fields.Nested(BaseAnnotationSchema, many=True)

# RNA Objects


# Protein Objects

class PolypeptideSchema(BaseRepositoryItemSchema):
    """
    Represents a protein
    """


class NucleaseSchema(PolypeptideSchema):
    """
    A protein that cuts DNA
    """


class RnaGuidedNucleaseSchema(NucleaseSchema):
    """
    The actual RNA Guided Nuclease itself.
    """


class RepositorySchema(Schema):
    """
    Represents a collection of biological objects either in vitro, such as
    an Inventory or in vivo such as a Genome
    """
    name = fields.String(required=True)