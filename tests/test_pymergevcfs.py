import unittest
from mergevcf.locations import *
from mergevcf.variantdict import *

class TestLocations(unittest.TestCase):

    def setUp(self):
        self.l1 = location('X',31,'+')
        self.l2 = location('Y',131,-1)
        self.l3 = location('Y',121,+1)
        self.l4 = location('Y',161,'-')
        self.ld = locationdict(window=40)
        self.ld[self.l1] = 'foo'
        self.ld[self.l2] = str(self.l2)
        self.ld[self.l3] = 'baz'

    def test_comparisons(self):
        self.assertTrue( self.l3 < self.l2 )
        self.assertTrue( self.l4 > self.l2 )
        self.assertFalse( self.l3.overlap(self.l4,strandTest=True,window=40) )
        self.assertTrue( self.l2.overlap(self.l4,strandTest=True,window=40) )

    def test_locationdict(self):
        self.assertTrue( self.l4 in self.ld )  # should overlap with l2
        self.assertTrue( self.ld.__str__() == "{(Y,121): 'baz', (Y,131,-L): '(Y,131,-L)', (X,31): 'foo'}")
        self.assertTrue( self.ld[self.l4].__str__() == "(Y,131,-L)")
        self.assertTrue( self.l3 in self.ld )
        self.assertTrue( self.l2 in self.ld )
        self.assertTrue( self.l1 in self.ld )

class TestVariantMap(unittest.TestCase):

    def setUp(self):
        self.l1 = location('X',31,'+')
        self.l2 = location('Y',131,-1)
        self.l3 = location('Y',121,+1)
        self.l4 = location('Y',161,'-')
        self.vmap = variantmap(0,40)
        self.vmap[(self.l1,('A','G'))] = 'sga'
        self.vmap[(self.l2,('ACGTA','A'))] = 'sanger'
        self.vmap[(self.l1,('A','G'))] = 'dkfz'
        self.vmap[(self.l2,self.l3)] = 'someone_else'
        self.vmap[(self.l4,self.l3)] = 'different'
        self.vmap[(self.l2,None)] = 'loosend'

    def test_vardict(self):
        self.assertTrue( (self.l1,('A','G')) in self.vmap )
        self.assertTrue( (self.l1,('C','G')) not in self.vmap )
        self.assertTrue( (self.l2,self.l3) in self.vmap )

    def test_iteration(self):
        nin = 0
        for var in self.vmap:
            nin += 1
        self.assertTrue( nin == 4 )

if __name__ == '__main__':
    unittest.main()
