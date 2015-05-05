"""Exonerate Parser"""
import re
import logging

def parse_extended_vulgar(aln_file):
    """
    Parses one line extended vulgar alignment summaries 
    from exonerate output. Outputs alignment results 
    using in-house conventions for coordinates/strand
    etc.

    The extended vulgar (verbose useful labelled gapped 
    alignment report)[1] is a standard output format. 
    The extension is to include total-equivalenced bases 
    and mismatches.

    Example as it appears in output:
    extended-vulgar: query 0 36 + 1 39811 39847 + 180 36 0 M 36 36 

    This is specified to to exonerate using the ryo 
    (roll your own) input:
    ryo='extended-vulgar: %S %et %em %V'

    The %et and %em expand to equivalenced bases and 
    mismatches. The '%S' expands to standard sugar 
    format, which is

    query_id
    query_start
    query_end
    query_strand
    target_id
    target_start
    target_end
    target_strand
    score

    The %V expands to a series of 3-tuples

       label query_length target_length 

    where label is one of

    M   Match (in the sense of equivalenced position, nt may differ)
    C   Codon 
    G   Gap 
    N   Non-equivalenced region 
    5   5' splice site 
    3   3' splice site 
    I   Intron 
    S   Split codon 
    F   Frameshift 

    The coordinates are on the forward strand and are 
    half open [,) zero based e.g. the substring 'el' 
    in 'hello' corresponds to [1,3). This follows the 
    in-house system so no shifts are needed. The 
    strand is given as '+' or '-' so this needs to be 
    converted to 1, -1. Coordinates on the reverse 
    strand are inverted, so have to be rectified for 
    output.

    Reference:
    [1] http://www.ebi.ac.uk/~guy/exonerate/exonerate.man.htm
    """

    # Regex
    line_re = r"^extended-vulgar:\s+(.*?)\s+(\d+)\s+(\d+)\s+([+-])\s+(.*?)\s+(\d+)\s+(\d+)\s+([+-])\s+(\d+)\s+(\d+)\s+(\d+)\s+(.*?)$"
    aln_re = r"([MCGN53ISF]+)\s+(\d+)\s+(\d+)"

    aln_fh = open(aln_file, 'r')
    aln_list = []

    for line in aln_fh:
        stripped = line.rstrip()
        if not line.startswith('extended-vulgar:'):
            continue
        match = re.match(line_re, stripped)
        if not match:
            raise Exception('Regex failed for ' + line)            
        aln_dict = {
            'query': match.group(1),
            'query_start': int(match.group(2)),
            'query_end': int(match.group(3)),
            'query_strand': 1 if match.group(4) == '+' else -1,
            'query': match.group(5),
            'target_start': int(match.group(6)),
            'target_end': int(match.group(7)),
            'target_strand': 1 if match.group(8) == '+' else -1,
            'score': int(match.group(9)),
            'equivalenced': int(match.group(10)),
            'mismatches': int(match.group(11)),
            'detail': []
        }
        alignment_detail = match.group(12)
        # Rectify coordinates if -ve strand
        for prefix in 'query', 'target':
            if aln_dict[prefix + '_strand'] == -1:
                ( aln_dict[prefix + '_start'], aln_dict[prefix + '_end'] ) = \
                ( aln_dict[prefix + '_end'], aln_dict[prefix + '_start'] )
        if not aln_dict[prefix + '_end'] > aln_dict[prefix + '_start']:
            raise Exception('Invalid interval ' + str(aln_dict))
        # Switch directions to get the query +1
        if aln_dict['query_strand'] == -1:
            aln_dict['target_strand'] *= -1
            aln_dict['query_strand'] = 1
        # Unpack the alignment
        match = re.findall(aln_re, alignment_detail)
        if not len(match):
            raise Exception('Regex failed for alignment ' + alignment_detail)            
        for triplet in match:
            triplet_dict = dict(zip(
                ('label', 'query_length', 'target_length'), 
                triplet))
            aln_dict['detail'].append(triplet_dict)

        # Append!
        aln_list.append(aln_dict)

    return aln_list


if __name__ == "__main__":
    pass






