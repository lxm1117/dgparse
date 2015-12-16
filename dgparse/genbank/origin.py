import string
from dgparse import sequtils

def recurse_origin(line, lines, out):
    columns = line.split()
    if not columns:
        return recurse_origin(next(lines), lines, out)
    bases = ''.join(map(string.upper, columns[1:]))
    if columns[0] == '//' or (not columns[0].isdigit() and not set(bases).issubset(sequtils.AMBIG_CHAR)):
        return line, lines, out
    else:
        out['origin'] = ''.join([out['origin'], bases])
        return recurse_origin(next(lines), lines, out)
