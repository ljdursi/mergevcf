"""
Define a variant map - a set of simple substitutions (SNVs, INDELs)
in one LocationPairDict, and a set of SVs in another
"""
from locations import Location, LocationDict
import vcf
import vcftobreakpoints as svvcf

def __checkvalidpairlocs__(candidate):
    if type(candidate) is not tuple:
        return False
    if len(candidate) != 2:
        return False
    if type(candidate[0]) is Location and type(candidate[1]) is Location:
        return True
    return False

class LocationPairDict(object):
    """
    Builds on locations.LocationDict, creates a dict of dicts
    for a given Location _pairs_
        dict[Location1][Location2] -> list of entries
    """
    def __init__(self, window):
        self.__window = window
        self.__lpdict = LocationDict(self.__window)

    def __contains__(self, lpair):
        if not __checkvalidpairlocs__(lpair):
            raise KeyError("Required: tuple of locations")
        locn1 = lpair[0]
        locn2 = lpair[1]
        if not locn1 in self.__lpdict:
            return False
        return locn2 in self.__lpdict[locn1]

    def __getitem__(self, lpair):
        if type(lpair) is Location:
            return self.__lpdict[lpair]

        if not __checkvalidpairlocs__(lpair):
            raise KeyError("Required: tuple of locations")
        if not self.__contains__(lpair):
            raise KeyError("Not in location pair dicti")

        locn1 = lpair[0]; locn2 = lpair[1]
        return self.__lpdict[locn1][locn2]

    def __setitem__(self, lpair, entry):
        if not __checkvalidpairlocs__(lpair):
            raise KeyError("Required: tuple of locations")
        locn1 = lpair[0]
        locn2 = lpair[1]
        if not locn1 in self.__lpdict:
            self.__lpdict[locn1] = LocationDict(self.__window)
        if not locn2 in self.__lpdict[locn1]:
            self.__lpdict[locn1][locn2] = []
        self.__lpdict[locn1][locn2].append(entry)

    def keys(self):
        return self.__lpdict.keys()

    def __iter__(self):
        return self.__lpdict.__iter__()

#    def __next__(self):
#        return self.__lpdict.__next()


class VariantMap(object):
    """
    Define a variant map - a set of simple substitutions (SNVs, INDELs)
    in one LocationPairDict, and a set of SVs in another
    """
    def __init__(self, awindow, svwindow):
        self.__awindow = awindow
        self.__svwindow = svwindow

        # map locn -> allele (ref/alt)
        self.__alleledict = LocationDict(awindow)

        # map locn -> locn (for SVs - paired breakpoints)
        self.__svdict = LocationPairDict(svwindow)

        # map locn -> all the locations found which map to this one
        self.__locn1pos = LocationPairDict(svwindow)
        self.__locn2pos = LocationPairDict(svwindow)

        # all the records involving this location
        self.__records = LocationPairDict(svwindow)

    def __medianpos__(self, locn1, locn2):
        def median(l):
            mid = len(l) // 2
            ls = sorted(l)
            if len(ls) % 2 == 0:
                return (ls[mid-1] + ls[mid]) // 2
            else:
                return ls[mid]

        if not (locn1, locn2) in self.__svdict:
            return None, None
        posn1s = self.__locn1pos[(locn1, locn2)]
        posn2s = self.__locn2pos[(locn1, locn2)]
        return median(posn1s), median(posn2s)

    def __svpresent__(self, locn1, locn2):
        return (locn1, locn2) in self.__svdict

    def __allelepresent__(self, locn, allele):
        if not locn in self.__alleledict:
            return False
        return allele in self.__alleledict[locn]

    def __addsvcaller__(self, locn1, locn2, caller):
        self.__svdict[(locn1, locn2)] = caller

    def __addallelecaller__(self, locn, allele, caller):
        if not locn in self.__alleledict:
            self.__alleledict[locn] = {}
        if not allele in self.__alleledict[locn]:
            self.__alleledict[locn][allele] = []
        if not caller in self.__alleledict[locn][allele]:
            self.__alleledict[locn][allele].append(caller)

    def __contains__(self, vartuple):
        assert type(vartuple) is tuple or type(vartuple) is list
        assert len(vartuple) > 1
        locn = vartuple[0]
        other = vartuple[1]
        if type(other) is Location:
            return self.__svpresent__(locn, other)
        else:
            return self.__allelepresent__(locn, other)

    def __setitem__(self, vartuple, caller, record=None):
        """Add item to the appropriate dictionary, and update locations
        dictionary as well."""
        assert type(vartuple) is tuple or type(vartuple) is list
        assert len(vartuple) > 1
        locn = vartuple[0]
        other = vartuple[1]
        if other is None:
            other = Location(None, 0)
        if type(other) is Location:
            self.__addsvcaller__(locn, other, caller)
            self.__locn1pos[(locn, other)] = locn.__pos__
            self.__locn2pos[(locn, other)] = other.__pos__
            self.__records[(locn, other)] = (caller, record)
        else:
            self.__addallelecaller__(locn, other, caller)

    def __getitem__(self, vartuple):
        assert type(vartuple) is tuple or type(vartuple) is list
        assert len(vartuple) > 1
        if not self.__contains__(vartuple):
            raise KeyError(str(vartuple))
        else:
            locn = vartuple[0]
            other = vartuple[1]
            if type(other) is Location:
                return self.__svdict[locn][other]
            else:
                return self.__alleledict[locn][other]

    def __str__(self):
        output = ""
        for loc in self.__alleledict:
            chrom = loc.__chrom__
            pos = loc.__pos__
            for item in self.__alleledict[loc]:
                ref, alt = item
                output += "\t".join([chrom, str(pos), '.', ref, alt])
                output += "\tCallers="+",".join(self.__alleledict[loc][item])
                output += "\n"

        for loc1 in self.__svdict:
            chrom1 = loc1.__chrom__
            for loc2 in self.__svdict[loc1]:
                chrom2 = loc2.__chrom__
                pos1, pos2 = self.__medianpos__(loc1, loc2)
                output += "\t".join([chrom1, str(pos1), '.', 'n',
                                     'n['+chrom2+":"+str(pos2)+'['])
                output += "\tCallers="+",".join(self.__svdict[loc1][loc2])+"\n"
        return output

    def __repr__(self):
        return self.__str__()

    # forceSV is here to allow forcing a call
    # that looks like a huge indel to be treated as an SV
    def addrecord(self, record, caller="NA", forceSV=False):
        assert type(record) is vcf.model._Record

        if forceSV or (record.ALT is not None and len(record.ALT) > 0
          and type(record.ALT[0]) in [vcf.model._SV, vcf.model._Breakend]):
            bkptpairs = svvcf.breakpoints_fm_record(record)
            for pair in bkptpairs:
                self.__setitem__(pair, caller, record)
        else:
            for alt in record.ALT:
                if alt is None:
                    continue
                loc = Location(record.CHROM, int(record.POS))
                allele = (record.REF, str(alt))
                self.__setitem__((loc, allele), caller)

    def __iter__(self):
        def gen_iterator():
            for loc in self.__alleledict:
                for allele in self.__alleledict[loc]:
                    yield loc, allele, self.__alleledict[loc][allele]
            for loc1 in self.__svdict:
                for loc2 in self.__svdict[loc1]:
                    pos1, pos2 = self.__medianpos__(loc1, loc2)
                    yield loc1, loc2, self.__svdict[loc1][loc2], pos1, pos2, self.__records[(loc1, loc2)]
            raise StopIteration()

        return gen_iterator()
