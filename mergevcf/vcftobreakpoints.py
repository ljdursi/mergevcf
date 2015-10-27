#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re
import vcf
import locations as loc

__symbolicRE__ = re.compile(r'.*<([A-Z:]+)>.*')
__bpRE__ = re.compile(r'([ACGTNactgn\.]*)([\[\]])([a-zA-Z0-9\.]+:\d+)([\[\]])([ACGTNacgtn\.]*)')
__looseendRE__ = re.compile(r'[lL][oO][oO][sS][eE][eE][nN][dD]')

def stdchrom(chrom):
    """
    Strip 'chr' from chromosome if present.
    """
    if chrom[0] == 'c':
        return chrom[3:]
    else:
        return chrom

def order_breakpoints(loc1, loc2):
    """
    Normalizes two breakpoints by ordering them lexicographically,
    flipping the half-interval between them if necessary, 
    and ensuring the first is on the positive strand.
    """
    # already ordered?  Just return them
    if loc1 < loc2:
        return loc1, loc2

    # Otherwise, swap and either flip strands (if different) or extents
    first, second = loc2, loc1
    if first.isRC():
        first, second = first.rc(), second.rc()

    return (first, second)

def otherPosnSymbolic(info):

    def getIfPresent(name, infodict):
        res = None
        if name in infodict:
            res = infodict[name]
            if type(res) is list:
                res = res[0]
        return res

    fields = ['CHR2', 'END', 'CT', 'SVTYPE', 'SVLEN']
    results = []
    for field in fields:
        results.append(getIfPresent(field, info))

    # handle chr2 specially
    if results[0] is not None:
        results[0] = stdchrom(results[0])

    # svclass overrides svtype
    svclass = getIfPresent('SVCLASS', info)
    if svclass is not None:
        results[3] = svclass

    transdict = {"inversion":"INV", "deletion":"DEL",
            "tandem_dup":"DUP:TANDEM", "duplication":"DUP",
            "long_range":"TRA", "inter_chr":"TRA"}

    if results[3] in transdict:
        results[3] = transdict[results[3]]

    return tuple(results)

def ctAndLocFromBkpt(ref, pre, delim1, pair, delim2, post):
    """
    extract strand/orientation, position, and (possibly) inserted string
    length from record eg [9:1678896[N.
    result is result from re.match with the explicit BP regexp
    """

    chpos = pair.split(':')
    chr2 = stdchrom(chpos[0])
    pos2 = int(chpos[1])
    assert delim1 == delim2  # '['..'[' or ']'...']'
    joined_after = True
    extend_right = True
    connect_sequence = ""

    if len(pre) > 0:
        connect_sequence = pre
        joined_after = True
        assert len(post) == 0
    elif len(post) > 0:
        connect_sequence = post
        joined_after = False

    if delim1 == "]":
        extend_right = False
    else:
        extend_right = True

    indellen = len(connect_sequence) - len(ref)

    if joined_after:
        if extend_right:
            connection_type = '3to5'
        else:
            connection_type = '3to3'
    else:
        if extend_right:
            connection_type = '5to5'
        else:
            connection_type = '5to3'

    return connection_type, chr2, pos2, indellen

def translocation(chr1, pos1, chr2, pos2, connection_type):
    """
    Returns a pair of locations corresponding to the two chr/pos pairs
    and the connection type.  Deterimines the strands and half intervals
    based on the connection type (one of '3to5','5to3', '3to3', '5to5')
    which is output by many callers.
    """
    second_strand = "+"
    first_extend_right = False
    second_extend_right = True
    if connection_type is None or connection_type == '3to5':
        pass
    elif connection_type == '5to3':
        first_extend_right, second_extend_right = (True, False)
    elif connection_type == '5to5':
        first_extend_right, second_extend_right, second_strand = (True, True, '-')
    elif connection_type == '3to3':
        first_extend_right, second_extend_right, second_strand = (False, False, '-')

    return (loc.Location(chr1, pos1, "+", first_extend_right),
            loc.Location(chr2, pos2, second_strand, second_extend_right))

#TODO: refactor
def breakpoints_fm_record(record):
    """Returns a list of pair(s) of breakpoints corresponding to the record."""
    chr1, pos1 = stdchrom(record.CHROM), int(record.POS)
    first = loc.Location(chr1, pos1)

    if record.ALT is None:
        return [(first, loc.Location(None, 0, "+", True))]

    ref = str(record.REF)
    breakpoint_pairs = []

    for alt in record.ALT:
        if alt is None:
            altstr = "."
        else:
            altstr = str(alt)

        # get all available information from the record
        chr2, pos2, connection_type, svtype, svlen = otherPosnSymbolic(record.INFO)

        # defaults
        if chr2 is None:
            chr2 = chr1

        # look for symbolic SVTYPE information in the alt field (eg, <DEL>)
        resultSym = re.match(__symbolicRE__, altstr)
        if resultSym:
            svtype = resultSym.group(1)

        # explicit BP - look for chr2 and CT information from the alt field
        # (eg, N[chr2:123123[)
        resultBP = re.match(__bpRE__, altstr)
        if resultBP:
            connection_type, chr2, pos2, indellen = ctAndLocFromBkpt(str(record.REF),
                                           resultBP.group(1), resultBP.group(2),
                                           resultBP.group(3), resultBP.group(4),
                                           resultBP.group(5))
            if svlen is None:
                svlen = indellen

        # looseend; no paired BP
        resultLE = __looseendRE__.search(str(record.FILTER))
        if altstr == ref+"." or resultLE:
            chr2, pos2 = (None, 0)
            if connection_type is None:
                connection_type = '5to3'
        if altstr == "."+ref:
            chr2, pos2 = (None, 0)
            if connection_type is None:
                connection_type = '3to5'

        # if nothing else, treat as an indel
        if not svtype and not resultBP and not resultLE and not resultSym:
            reflen = len(ref)
            if pos2 is None:
                pos2 = pos1 + reflen
            if svlen is None:
                svlen = len(altstr) - reflen
            if svtype is None:
                if svlen > 0:
                    svtype = 'INS'
                else:
                    svtype = 'DEL'

        if svtype == "INV":
            assert chr1 == chr2
            breakpoint_pairs.append(translocation(chr1, pos1,   chr1, pos2-1, '3to3'))
            breakpoint_pairs.append(translocation(chr1, pos1+1, chr1, pos2,   '5to5'))
        elif svtype in ["DUP:TANDEM", "DUP"]:
            if chr2 is None:
                chr2 = chr1
            if connection_type is None or connection_type == "5to3" or connection_type == "3to5":
                if pos1 > pos2:
                    pos1, pos2 = pos2, pos1
                breakpoint_pairs.append(translocation(chr1, pos1, chr2, pos2, "5to3"))
            else:
                breakpoint_pairs.append(translocation(chr1, pos1, chr2, pos2, connection_type))
        elif svtype in ["INS", "INS:ME:L1"]:
            assert chr1 == chr2
            breakpoint_pairs.append((first, loc.Location(chr1, pos1+1, "+", True)))
        elif svtype in ["DEL:ME:ALU", "DEL"]:
            assert chr1 == chr2
            breakpoint_pairs.append(translocation(chr1, min(pos1, pos2), chr1,
                                           max(pos1, pos2), '3to5'))
        else:
            # translocation is the default
            if connection_type is None:
                connection_type = "3to5"
            if (svtype is not None) and (not svtype in ["TRA", "BND"]):
                print("Got unknown record of type",
                        svtype, altstr, str(record), file=sys.stderr)
                print("Hoping for best", file=sys.stderr)
            breakpoint_pairs.append(translocation(chr1, pos1, chr2, pos2, connection_type))

    return [order_breakpoints(bp1, bp2) for bp1, bp2 in breakpoint_pairs]

def add_bkpt_to_dict_dict(bkpt, ld):
    """ Add to dictionary if not none """
    if not bkpt is None:
        ld[bkpt] = 1

def breakpoint_in_dict_dict(bkpt, ldict):
    """ Lookup in dictionary if not None """
    if not bkpt is None:
        return bkpt in ldict
    else:
        return False

def vcftobkpts(infile, outfile, width):
    """
    Routine which parses a VCF file and outputs a minimal list of breakpoints,
    for testing
    """
    firstbkpts = loc.LocationDict(width)
    pairbkpts = loc.LocationDict(width)

    reader = vcf.Reader(infile)
    for record in reader:
        if record.FILTER == "PASS" or record.FILTER == "." or record.FILTER is None or (type(record.FILTER) is list and len(record.FILTER) == 0):
            breakpoint_pairs = breakpoints_fm_record(record)
            for pair in breakpoint_pairs:
                add_bkpt_to_dict_dict(pair[0], firstbkpts)
                add_bkpt_to_dict_dict(pair[1], pairbkpts)

    # count how many breakpoints weren't in both dicts, for diagnostics;
    # then add everything to first dict for outputting
    nunmatched = 0
    for bkpt in firstbkpts:
        if not breakpoint_in_dict_dict(bkpt, pairbkpts):
            nunmatched += 1

    for bkpt in pairbkpts:
        if not breakpoint_in_dict_dict(bkpt, firstbkpts):
            nunmatched += 1
            add_bkpt_to_dict_dict(bkpt, firstbkpts)

    print("#Num breakpoints not in both lists:", nunmatched, file=sys.stderr)

    # now output everything
    for bkpt in firstbkpts:
        chrom, pos, _, _ = bkpt.asTuple()
        start = pos-width/2
        if start < 0:
            start = 0
        print("{0}	{1}	{2}".format(chrom, start, pos+width/2), file=outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
                         default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
                         default=sys.stdout)
    defaultwidth = 300
    parser.add_argument('-w', '--width', type=int,
                        help="width of breakpoint region: default("+str(defaultwidth)+")",
                        default=defaultwidth)

    args = parser.parse_args()
    sys.exit(vcftobkpts(args.infile, args.outfile, args.width))
