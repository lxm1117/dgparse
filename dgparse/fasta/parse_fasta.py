# encoding=utf-8
"""
Parse a Fasta file.
"""
import hashlib
import logging


from exceptions import UnicodeDecodeError
from exceptions import TypeError
from dgparse.sequtils import DNA_CHAR
from dgparse.sequtils import NOT_DNA

from dgparse.exc import ParserException

log = logging.getLogger(__name__)


def parse_fasta_str(fasta_str, default={}):
    """
    Parse a fasta formatted String (not a file)
    :param fasta_str:
    :param default:
    :return:
    """
    if not isinstance(fasta_str, basestring):
        raise TypeError()
    lines = fasta_str.split("\n")
    records = {}
    result = []
    count = 0
    current_header = None
    for line in lines:
        if "rtf" in line:
            msg = "This string is RTF formatted"
            log.warn(msg)
            raise ParserException(msg)
        if line and ">" in line:
            seqrec = dict(default)  # make new seqrec
            try:
                seqrec['description'] = line
            except UnicodeDecodeError:
                message = "Unicode Decode Error on line: {}".format(line)
                log.warn(message)
                raise ParserException(message)
            seqrec['file_format'] = u'fasta'
            seqrec['sequence'] = []
            seqrec['is_circular'] = False if "linear" in line else True
            try:
                name = parse_fasta_header_line(line)
            except IndexError:
                assert not line, repr(line)
                # Should we use SeqRecord default for no ID?
                name = u"Plasmid"
                if count > 0:
                    name = name + "-" + unicode(count)
                count += 1
            seqrec['name'] = name
            records[line] = seqrec  # index by header line
            current_header = line
        elif line and (line[0] in DNA_CHAR):
            clean_line = line.strip('{}/ \n\\')
            if not current_header:
                message = ('This file does not containa valid FASTA header '
                           'line. Please include >Name on the first line of '
                           'the file or upload this sequence as Plain Text.')
                raise ParserException(message)

            records[current_header]['sequence'].append(clean_line)

        else:
            pass
    # Now process the header-body pairs into sequence records
    for seqrec in records.itervalues():
        lines = seqrec['sequence']
        seqstr = unicode("".join(lines).replace(" ", "").replace("\r",
                                                                 "").upper())
        if NOT_DNA.search(seqstr):
            raise ParserException('Invalid sequence.')
        new_sequence = {
            'sha1': hashlib.sha1(seqstr).hexdigest(),
            'seq': seqstr,
        }
        seqrec['sequence'] = new_sequence
        seqrec['length'] = len(seqstr)
        result.append(seqrec)
    if not result:
        raise ParserException('Nothing could be parsed')
    return result


def parse_fasta_file(fasta_file, default={}):
    """
    Parse a fasta file uploaded rather than a fasta formatted string pasted into
    a form field.
    The result of this should be a de-duplicated list of dictionaries
    """
    records = {}
    result = []
    count = 0
    current_header = None
    saved_contents = list()
    for line in fasta_file.readlines():
        saved_contents.append(line)
        if "rtf" in line:
            return None
        if line == "":
            return  # Premature end of file, or just empty?
        if line and ">" in line:
            seqrec = dict(default)  # make new seqrec
            seqrec['description'] = unicode(line, 'utf-8')
            seqrec['file_format'] = u'fasta'
            seqrec['sequence'] = []
            if "linear" in line:
                seqrec['is_circular'] = False
            if 'circular' in line:
                seqrec['is_circular'] = True
            try:
                name = parse_fasta_header_line(line)
            except IndexError:
                assert not line, repr(line)
                name = fasta_file.name.split("/")[-1]
                if count > 0:
                    name = name + "-" + str(count)
                count += 1
            seqrec['name'] = unicode(name, 'utf-8')
            records[line] = seqrec  # index by header line
            current_header = line
        elif line and (line[0] in DNA_CHAR):

            if "N" in line:
                log.warn("Ambiguous Base found in this sequence. Removing")
            if not current_header:
                msg = "This file does not containa valid FASTA header line.\
                       Please include >Name on the first line of the file."
                raise ParserException(msg)

            clean_line = line.translate(None, '{}/ \n\\')
            if not current_header:
                raise ParserException("The Current Header is NONE")

            records[current_header]['sequence'].append(clean_line)
            # Remove trailing whitespace, and any internal spaces
            # (and any embedded \r which are possible in mangled files
            # when not opened in universal read lines mode)
        else:
            pass
    if len(records) == 1:
        records[current_header]['file_contents'] = '\n'.join(saved_contents)
    # Now process the header-body pairs into sequence records
    for seqrec in records.itervalues():
        lines = seqrec['sequence']
        seqstr = unicode("".join(lines).replace(" ", "").replace("\r",
                                                                 "").upper())
        if NOT_DNA.search(seqstr):
            msg = u"Invalid character found in {n}".format(n=seqrec[u'name'])
            raise ParserException(msg)
        else:
            new_sequence = {
                'seq': seqstr,
                'sha1': hashlib.sha1(seqstr).hexdigest(),
            }
        seqrec['sequence'] = new_sequence
        seqrec['length'] = len(seqstr)
        result.append(seqrec)
    return result


def parse_fasta_header_line(line):
    """
    Parse the header line.
    :param line:
    :return:
    """
    name = line.split(" ")
    #print "Split results:", name
    return name[0][1:]
