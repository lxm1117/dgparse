
import hashlib
import logging


from dgparse.utility import DNA_CHAR
from dgparse.utility import NOT_DNA
from dgparse.exc import ParserException

log = logging.getLogger(__name__)


def parse_fasta_header_line(line):

    name = line.split(" ")
    #print "Split results:", name
    return name[0][1:]


def parse_fasta_file(fasta_file, default={}):
    """
    Parse a fasta file uploaded rather than a fasta formatted string pasted into
    a form field.
    The result of this should be a de-duplicated list of dictionaries

    """
    records = {}
    result = []
    count = 0
    #print "This file name is:", fasta_file.name
    current_header = None
    saved_contents = list()
    for line in fasta_file.readlines():
        saved_contents.append(line)
        if "rtf" in line:
            print "{} is RTF formatted".format(fasta_file.name)
            return None
        if line == "" :
            return #Premature end of file, or just empty?
        if line and ">" in line:
            #print "We have a new record:", line
            seqrec = dict(default) #make new seqrec
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
                #Should we use SeqRecord default for no ID?
                name = fasta_file.name.split("/")[-1]
                if count > 0:
                    name = name + "-" + str(count)
                count += 1
            #print "Name:", name
            seqrec['name'] = unicode(name, 'utf-8')
            records[line] = seqrec #index by header line
            current_header = line
        elif line and (line[0] in DNA_CHAR):

            if "N" in line:
                log.warn("Ambiguous Base found in this sequence. Removing")
            clean_line = line.translate(None, '{}/ N\n\\')
            #print "Clean Line:", clean_line
            if not current_header:
                msg = "This file does not containa valid FASTA header line.\
                       Please include >Name on the first line of the file."
                raise ParserException(msg)
                #raise Exception("The Current Header is NONE")

            clean_line = line.translate(None, '{}/ \n\\')
            #print "Clean Line:", clean_line
            if not current_header:
                raise ParserException("The Current Header is NONE")

            records[current_header]['sequence'].append(clean_line)
            #records[name].append(line.rstrip('}{\\/'))
            #Remove trailing whitespace, and any internal spaces
            #(and any embedded \r which are possible in mangled files
            #when not opened in universal read lines mode)
        else:
            #print "Dirty line:", line
            pass
    if len(records) == 1:
        records[current_header]['file_contents'] = '\n'.join(saved_contents)
    #Now process the header-body pairs into sequence records
    for seqrec in records.itervalues():
        lines = seqrec['sequence']
        seqstr = unicode("".join(lines).replace(" ", "").replace("\r",
                                                                 "").upper())
        if NOT_DNA.search(seqstr):
            new_sequence = None
            msg = u"Invalid character found in {n}".format(n=seqrec[u'name'])
            raise ParserException(msg)
        else:
            new_sequence = {
                'bases': seqstr,
                'sha1': hashlib.sha1(seqstr).hexdigest(),
            }
        seqrec['sequence'] = new_sequence
        seqrec['length'] = len(seqstr)
        result.append(seqrec)
    return result



