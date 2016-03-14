"""
Classes for signed genomic locations - eg, endpoints of half-intervals
"""

class location(object):
    __strands = {+1:True, -1:False, "+":True, "-":False, True:True, False:False}
    __strandstrs = {True:"+", False:"-"}
    __extendstrs = {True:"R", False:"L"}
    def __init__(self, chrom, pos, strand=True, extendsRight=False):
        assert type(pos) is int
        assert strand in self.__strands
        assert type(extendsRight) is bool
        self.__chrom__ = str(chrom)
        self.__pos__ = pos
        self.__strand__ = self.__strands[strand]
        self.__right__ = extendsRight

    def __hash__(self):
        return hash(self.__chrom__) ^ hash(self.__pos__) ^ hash(self.__strand__) ^ hash(self.__right__)
        
    def __cmp__(self, other):
        """Compare location of two positions, regardless of strand or direction."""
        def _strsAsIntsIfPossible(str1, str2):
            try:
                i1 = int(str1)
                i2 = int(str2)
            except ValueError:
                return str1, str2
            return i1, i2

        if other is None:
            return -1 
        if self.__chrom__ != other.__chrom__:
            ic1, ic2 = _strsAsIntsIfPossible(self.__chrom__, other.__chrom__)
            return cmp(ic1, ic2)
        return cmp(self.__pos__, other.__pos__)

    def overlap(self, other, strandTest=False, window=0):
        """Compare location of two endpoints with the given window, optionally
        testing for strand"""
        if other is None:
            return False
        if self.__chrom__ != other.__chrom__:
            return False
        if strandTest and (self.__strand__ != other.__strand__):
            return False
        return abs(self.__pos__ - other.__pos__) <= window

    def __add__(self, offset):
        """Add offset to position"""
        return location(self.__chrom__, self.__pos__+offset, self.__strand__, self.__right__)

    def __str__(self):
        """Convert to string"""
        return '(%s,%d,%s%s)' % (self.__chrom__, self.__pos__, self.__strandstrs[self.__strand__], self.__extendstrs[self.__right__])

    def __repr__(self):
        """Representation for printing"""
        extrastr=""
        if (not self.__strand__) or (self.__right__):
            extrastr = ",%s%s" % (self.__strandstrs[self.__strand__], self.__extendstrs[self.__right__])
        return '(%s,%d%s)' % (self.__chrom__, self.__pos__, extrastr)

    def rc(self):
        """Reverse strand of current strand"""
        return location(self.__chrom__, self.__pos__, not self.__strand__, self.__right__)

    def isRC(self):
        """Is on reverse strand?"""
        return self.__strand__ == False

    def switchExtent(self):
        return location(self.__chrom__, self.__pos__, self.__strand__, not self.__right__)

    def asTuple(self):
        return self.__chrom__, self.__pos__, self.__strand__, self.__right__

    def withPos(self, pos):
        """Replace position of location"""
        return location(self.__chrom__, pos, self.__strand__, self.__right__)

    @property
    def pos(self):
        return self.__pos__

    @property
    def chrom(self):
        return self.__chrom__

class locationdict(dict):
    def __init__(self, window, *args, **kwargs):
        self.__window = window
        self.__search = [0] + [item for pm in zip(range(1,window+1),range(-1,-window-1,-1)) for item in pm]
        super(locationdict,self).__init__(*args, **kwargs)

    def keys(self):
        return super(locationdict, self).keys()

    def values(self):
        return [self[key] for key in self]

    def itervalues(self):
        return (self[key] for key in self)

    def __find__(self, locn):
        if not type(locn) is location:
            raise ValueError("Not Location: "+locn.__str__())
        present = False
        foundoff = None
        for off in self.__search:
            if super(locationdict, self).__contains__(locn+off):
                present = True
                foundoff = off
        return present, foundoff

    def __contains__(self, locn):
        if not type(locn) is location:
            raise ValueError("Not Location: "+locn.__str__())
        present, foundloc = self.__find__(locn)
        return present

    def __getitem__(self, locn):
        if not type(locn) is location:
            raise ValueError("Not Location: "+locn.__str__())
        present, foundoff = self.__find__(locn)
        if not present:
            raise KeyError(locn.__str__())
        else:
            return super(locationdict, self).__getitem__(locn + foundoff)
