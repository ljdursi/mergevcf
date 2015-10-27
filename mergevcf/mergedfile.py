"""
Functions to operate on merged VCF files and to perform the merge
of multiple VCF files
"""
import vcf
import variantdict

def is_mapped_to_chromosome(chrom):
    "Is the variant on a human chromosome or another contig?"
    if chrom[0:2] in ["GL", "MT", "hs"] or chrom[0:1] == "M":
        return False
    return True

def int_if_possible(stringval):
    "Returns int representation of string if available"
    if type(stringval) is list:
        stringval = stringval[0]
    try:
        i = int(stringval)
    except ValueError:
        i = stringval
    return i

def make_info_dict(records, pos1, pos2):
    """generate 'median' results for info fields from all records"""
    fields = ['CHR2', 'END', 'SVTYPE', 'SVLEN']
    info = {}
    for field in fields:
        answers = []
        for record in records:
            if field in record.INFO:
                answers.append(int_if_possible(record.INFO[field]))
        nanswers = len(answers)
        if nanswers > 0:
            sortedanswers = sorted(answers)
            medianpos = nanswers/2
            mediananswer = sortedanswers[medianpos]
            if not mediananswer == 0:
                info[field] = mediananswer

    if 'SVTYPE' in info and info['SVTYPE'] in ['DUP', 'DUP:TANDEM', 'DEL', 'INV']:
        if not 'SVLEN' in info:
            info['SVLEN'] = pos2-pos1
    return info

def bkpt_ref_alt_from_pair(loc1, loc2, refstr="N"):
    """Return the breakpoint-style ref-alt given two locations"""
    altafter = loc1.__right__ == False

    if loc2 is None or loc2.__chrom__ is None:
        bkptstr = "."
    else:
        if loc2.__right__:
            delim = "["
        else:
            delim = "]"

        assert loc2.__strand__ == (loc1.__right__ != loc2.__right__)
        bkptstr = "%s%s:%d%s" % (delim, loc2.__chrom__, loc2.__pos__, delim)

    if altafter:
        altstr = "%s%s" % (refstr, bkptstr)
    else:
        altstr = "%s%s" % (bkptstr, refstr)

    return refstr, altstr

#TODO: refactor
def merge(filenames, programs, forceSV, outfile, slop=0, verbose=True,
          filter_chrom=True, no_filter_variants=False):
    """Merge several VCFs from different programs into a new VCF file."""

    def passed_variant(record):
        """Returns true if the variant is PASS in the VCF file, 
           or if we are not filtering"""
        return record.FILTER is None or len(record.FILTER) == 0 or no_filter_variants

    def info_string(callers, infodict):
        """Generate an info string for merged records"""
        infostring = ""
        fields = ['CHR2', 'END', 'SVTYPE', 'SVLEN']
        for field in fields:
            if field in infodict:
                res = infodict[field]
                if type(res) is list:
                    res = res[0]
                infostring = infostring + ';'+field+'='+str(res)
        return "Callers="+",".join(list(set(callers)))+infostring

    calldict = variantdict.VariantMap(awindow=0, svwindow=slop)
    for (infile, program) in zip(filenames, programs):
        count = 0
        try:
            vcf_reader = vcf.Reader(open(infile, 'r'))
            for record in vcf_reader:

                if not passed_variant(record):
                    continue

                if filter_chrom and not is_mapped_to_chromosome(record.CHROM):
                    continue

                if verbose:
                    if count == 0:
                        print record, program
                    count += 1
                    if count == 100:
                        count = 0

                calldict.addrecord(record, program, forceSV)
        except (RuntimeError, TypeError, NameError, AttributeError):
            pass

    # Write the results in a master vcf file for the sample
    outfile.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for variant in calldict:
        if len(variant) == 3:   # snv/indel
            loc, allele, callers = variant
            if allele is None:
                print "Allele is none: loc, allele, callers = ", loc, allele, callers
                continue
            chrom, pos, _, _ = loc.asTuple()
            fields = [chrom, str(pos), allele[0], allele[1]]
            vcfline = "\t".join([fields[0], fields[1], ".", fields[2],
                                 fields[3], "255", "PASS",
                                 "Callers=" + ",".join(callers)])
            outfile.write(vcfline + "\n")
        else:
            loc1, loc2, callers, medianpos1, medianpos2, recordscalled = variant
            if filter_chrom and not is_mapped_to_chromosome(loc2.chrom):
                continue

            records = [r for _, r in recordscalled]

            avgloc1 = loc1.withPos(medianpos1)
            avgloc2 = loc2.withPos(medianpos2)
            ref, alt = bkpt_ref_alt_from_pair(avgloc1, avgloc2)
            vcfline = "\t".join([avgloc1.__chrom__, str(avgloc1.__pos__), '.',
                                 ref, alt, '255', 'PASS',
                                 info_string(callers, make_info_dict(records, medianpos1, medianpos2))])
            outfile.write(vcfline + "\n")
            for caller, rec in recordscalled:
                outfile.write("#"+str(rec)+" ("+caller+")\n")

    outfile.close()

def read_merged_calls(infile, filter_chrom=True, skipcallers=None):
    """Read a merged callset, and return:
        - dictionary: caller name -> caller idx
        - callsets(list of lists): [calleridx][callidx]
        - calls: callidx -> record from merged"""
    invcf = vcf.Reader(infile)
    callerIdx = 0
    callIdx = 0
    callsets = []
    callIdxToCall = []
    callerIdxDict = {}

    if skipcallers is None:
        skipcallers = []

    for rec in invcf:
        ncalledthis = 0
        if filter_chrom and not is_mapped_to_chromosome(rec.CHROM):
            continue
        callers = [c for c in rec.INFO['Callers'] if not c in skipcallers]
        called = []
        for caller in callers:
            if not (caller in called) and not (caller in skipcallers):
                called.append(caller)

                if not caller in callerIdxDict:
                    callerIdxDict[caller] = callerIdx
                    callerIdx += 1
                    callsets.append([])

                callsets[callerIdxDict[caller]].append(callIdx)
                ncalledthis += 1

        assert len(called) == ncalledthis
        if ncalledthis > 0:
            chrom = rec.CHROM
            posstart = rec.POS
            callIdxToCall.append((len(called), chrom, posstart,
                                  str(rec.REF), str(rec.ALT[0]),
                                  ",".join(called)))
            callIdx += 1

    return callerIdxDict, callsets, callIdxToCall

