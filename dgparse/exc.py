"""
Parser Exceptions shared by all Parsers
"""

class ParserException(Exception):
    pass

class FormatException(ParserException):
    """
    Use for badly formatted files
    """
    pass
