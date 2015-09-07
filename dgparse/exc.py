"""
Parser Exceptions shared by all Parsers
"""

from marshmallow import ValidationError


class ParserException(Exception):
    pass


class FormatException(ParserException):
    """
    Use for badly formatted files
    """
    pass


class UndefinedRecordType(FormatException):
    """
    Use when no record type is stated or inferred.
    """


class NonRecordEntry(FormatException):
    """
    Use when a CSV or Excel file contains extraneous records
    """


class IllegalCharacter(ValidationError):
    """
    Raise when a disallowed character is present
    """


class NoSequence(ValidationError):
    """
    Raise when an empty sequence or pattern is provided
    """


class NullCoordinates(ParserException):
    """Features must have a length of 1 or more"""
    pass
