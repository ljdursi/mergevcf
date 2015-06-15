#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re
import vcf
import locations as loc

__symbolicRE__ = None
__bpRE__ = None
__looseendRE__ = None

def setupREs():
    global __symbolicRE__
    global __bpRE__
    global __looseendRE__
    if __symbolicRE__ is None or __bpRE__ is None or __looseendRE__ is None:
        __symbolicRE__ = re.compile(r'.*<([A-Z:]+)>.*')
        __bpRE__ = re.compile(r'([ACGTNactgn\.]*)([\[\]])([a-zA-Z0-9\.]+:\d+)([\[\]])([ACGTNacgtn\.]*)')
        __looseendRE__ = re.compile(r'[lL][oO][oO][sS][eE][eE][nN][dD]')

def stdchrom(chrom):
    if chrom[0]=='c':
        return chrom[3:]
    else:
        return chrom

def orderBreakpoints(loc1, loc2):
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

    fields = ['CHR2','END','CT','SVTYPE','SVLEN']
    results = []
    for field in fields:
        results.append( getIfPresent(field, info) )

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
    # extract strand/orientation, position, and (possibly) inserted string
    # length from record eg [9:1678896[N.
    # result is result from re.match with the explicit BP regexp

    chpos=pair.split(':')
    chr2 = stdchrom(chpos[0])
    pos2 = int(chpos[1])
    assert delim1 == delim2  # '['..'[' or ']'...']'
    joinedAfter = True
    extendRight = True
    connectSeq = ""

    if len(pre) > 0:
        connectSeq = pre 
        joinedAfter = True
        assert len(post) == 0
    elif len(post) > 0:
        connectSeq = post
        joinedAfter = False

    if delim1 == "]":
        extendRight = False
    else:
        extendRight = True

    indellen = len(connectSeq) - len(ref)

    if joinedAfter:
        if extendRight:
            ct = '3to5'
        else:
            ct = '3to3'
    else:
        if extendRight:
            ct = '5to5'
        else:
            ct = '5to3'
    
    return ct, chr2, pos2, indellen

def translocation(chr1, pos1, chr2, pos2, ct):
    secondStrand = "+"
    firstExtendRight = False
    secondExtendRight = True
    if ct is None or ct == '3to5':
        pass
    elif ct == '5to3':
        firstExtendRight = True; secondExtendRight = False
    elif ct == '5to5':
        firstExtendRight = True; secondExtendRight = True; secondStrand = '-'
    elif ct == '3to3':
        firstExtendRight =False; secondExtendRight =False; secondStrand = '-'

    return ( loc.location(chr1, pos1, "+",          firstExtendRight), 
             loc.location(chr2, pos2, secondStrand, secondExtendRight) )

def breakpointsFromRecord(record):
    """Returns a list of pair(s) of breakpoints corresponding to the record."""
    global __symbolicRE__
    global __bpRE__
    global __looseendRE__

    chr1, pos1 = stdchrom(record.CHROM), int(record.POS)
    first = loc.location(chr1, pos1)

    if record.ALT is None:
        return [(first, loc.location(None,0,"+",True))]

    setupREs()
    ref = str(record.REF)
    bkptPairs = []

    for alt in record.ALT:
        if alt is None:
            altstr = "."
        else:
            altstr = str(alt)

        # get all available information from the record 
        chr2, pos2, ct, svtype, svlen = otherPosnSymbolic(record.INFO)

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
            ct, chr2, pos2, indellen = ctAndLocFromBkpt(str(record.REF), resultBP.group(1),
                    resultBP.group(2), resultBP.group(3), resultBP.group(4),
                    resultBP.group(5))
            if svlen is None:
                svlen = indellen

        # looseend; no paired BP
        resultLE = __looseendRE__.search(str(record.FILTER))
        if altstr == ref+"." or resultLE:
            chr2 = None; pos2 = 0
            if ct is None:
                ct = '5to3'
        if altstr == "."+ref:
            chr2 = None; pos2 = 0
            if ct is None:
                ct = '3to5'

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
            bkptPairs.append( translocation(chr1, pos1,   chr1, pos2-1, '3to3') )
            bkptPairs.append( translocation(chr1, pos1+1, chr1, pos2,   '5to5') )
        elif svtype in ["DUP:TANDEM","DUP"]:
            if chr2 is None:
                chr2 = chr1
            if ct is None or ct=="5to3" or ct=="3to5":
                if pos1 > pos2:
                    pos1, pos2 = pos2, pos1
                bkptPairs.append( translocation(chr1, pos1, chr2, pos2, "5to3") )
            else:
                bkptPairs.append( translocation(chr1, pos1, chr2, pos2, ct) )
        elif svtype in ["INS","INS:ME:L1"]:
            assert chr1 == chr2
            bkptPairs.append((first, loc.location(chr1, pos1+1, "+", True)))
        elif svtype in ["DEL:ME:ALU","DEL"]:
            assert chr1 == chr2
            bkptPairs.append( translocation(chr1, min(pos1,pos2), chr1, max(pos1,pos2), '3to5') )
        else:
            # translocation is the default
            if ct is None:
                ct = "3to5"
            if (svtype is not None) and (not svtype in ["TRA","BND"]):
                print("Got unknown record of type", svtype, altstr, str(record), file=sys.stderr)
                print("Hoping for best",file=sys.stderr)
            bkptPairs.append( translocation(chr1, pos1, chr2, pos2, ct) )

    orderedPairs = [orderBreakpoints(bp1,bp2) for bp1,bp2 in bkptPairs]
    return orderedPairs

def addBkptToDictDict(bkpt, ld):
    if not bkpt is None:
        ld[bkpt] = 1

def bkptInDictDict(bkpt, ld):
    if not bkpt is None:
        return bkpt in ld
    else:
        return False

def vcftobkpts(infile, outfile, width):
    firstbkpts = loc.locationdict(width)
    pairbkpts  = loc.locationdict(width)

    reader = vcf.Reader(infile)
    for record in reader:
        if record.FILTER == "PASS" or record.FILTER == "." or record.FILTER is None or (type(record.FILTER) is list and len(record.FILTER) == 0):
            bkptPairs = breakpointsFromRecord(record)
            for pair in bkptPairs:
                addBkptToDictDict(pair[0], firstbkpts)
                addBkptToDictDict(pair[1], pairbkpts)

    # count how many breakpoints weren't in both dicts, for diagnostics;
    # then add everything to first dict for outputting
    nunmatched = 0
    for bp in firstbkpts:
        if not bkptInDictDict(bp, pairbkpts):
            nunmatched += 1

    for bp in pairbkpts:
        if not bkptInDictDict(bp, firstbkpts):
            nunmatched += 1
            addBkptToDictDict(bp, firstbkpts)

    print("#Num breakpoints not in both lists:",nunmatched,file=sys.stderr)

    # now output everything
    for bp in firstbkpts:
        chrom, pos, strand, extendsRight = bp.asTuple()
        start = pos-width/2
        if start < 0:
            start = 0
        print("{0}	{1}	{2}".format(chrom, start, pos+width/2), file=outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'), default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    defwidth=300
    parser.add_argument('-w','--width', type=int, help="width of breakpoint region: default("+str(defwidth)+")",default=defwidth)

    args = parser.parse_args()
    sys.exit( vcftobkpts(args.infile,args.outfile,args.width ) )
