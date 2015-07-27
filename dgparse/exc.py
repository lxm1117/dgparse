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