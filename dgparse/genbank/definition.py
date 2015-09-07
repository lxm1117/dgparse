from .constants import GENBANK_HEADERS


def recurse_definition(line, lines, out):
    columns = line.split()
    if not columns:
        return recurse_definition(next(lines), lines, out)
    if columns[0] != 'DEFINITION':
        if columns[0] in GENBANK_HEADERS:
            return line, lines, out
        out['definition'] = ' '.join([out['definition'], ' '.join(columns)])
        return recurse_definition(next(lines), lines, out)
    if columns[0] == 'DEFINITION':
        out['definition'] = ' '.join(columns[1:])
        return recurse_definition(next(lines), lines, out)

