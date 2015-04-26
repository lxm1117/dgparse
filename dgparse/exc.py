class ParserException(Exception):
    pass


class NoParserException(ParserException):
    pass


class FeaturesStopIteration(ParserException):
    pass


class NullCoordinates(ParserException):
    """Features must have a length of 1 or more"""
    pass