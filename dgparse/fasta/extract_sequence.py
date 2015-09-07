EMPTY = 'EMPTY'


def detect_circular(input_file):
    input_file.seek(0)
    first_line = input_file.readline()
    if 'linear' in first_line:
        return False
    input_file.seek(0)
    return True


def guess_format_from_content(content):
    record_format = None
    if not len(content):
        return EMPTY
    if content[0] == ">":
        record_format = "fasta"
    if "LOCUS" in content:
        record_format = "genbank"
    return record_format
