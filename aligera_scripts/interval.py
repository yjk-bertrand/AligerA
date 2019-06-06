#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
#  Class Interval that is a data structure use in the interval search tree.
__version__ = "1.3"



# =======================================================================================
#               IMPORTS
# =======================================================================================

from math import floor


# =======================================================================================
#               CLASS
# =======================================================================================


class Interval:
    """
    Data structure use in the interval search tree. It relies on the interval data 
    structure found at:
    https://algs4.cs.princeton.edu/93intersection/IntervalST.java.html
    This data structure is used in the STEP2 of the AligerA pipeline tool to trim
    the flanking regions of an alignment.
    """

    def __init__(self, lo, hi):

        if lo is None or hi is None:
            raise ValueError(
                "Interval: None endpoint not allowed: ({0}, {1})".format(lo, hi)
            )

        if lo <0 or hi <= 0:
            raise ValueError(
                "Interval: Negative lower endpoint and null upper endpoint are \
not allowed: ({0}, {1})".format(
                    lo, hi
                )
            )
            
                  
        if lo >= hi:
            raise ValueError(
                "Interval: Lower endpoint must be higher than upper \
endpoint: ({0}, {1})".format(
                    lo, hi
                )
            )

        self.lo = lo
        self.hi = hi


    @classmethod
    def from_tuple(cls, tuple_lo_hi):
        if len(tuple_lo_hi) != 2:
            raise ValueError(
                "Node: tuple can contain only two values {} ".format(tuple_lo_hi)
            )
        lo = tuple_lo_hi[0]
        hi = tuple_lo_hi[1]
        interval = cls(lo, hi)
        return interval

    
    def compare(self, other):
        #  ascending order of lo endpoint, breaking ties by hi endpoint
        if floor(self.lo) == floor(other.lo) and floor(self.hi) == floor(other.hi):
            return 0
        elif floor(self.lo) < floor(other.lo):
            return -1
        elif floor(self.lo) > floor(other.lo):
            return 1
        elif self.hi < other.hi:
            return -1
        else:  #  self.hi > other.hi
            return 1


    def contains(self, other):
        if self.lo <= other.lo and other.hi <= self.hi:
            return True
        return False


    def intersects(self, other):
        if self.hi < other.lo:
            return False
        if other.hi < self.lo:
            return False
        return True


    def length(self):
        return self.hi - floor(self.lo)


    def getMax(self):
        return self.hi


    def getMin(self):
        return self.lo


    def __repr__(self):
        return "({0}, {1})".format(self.lo, self.hi)


    def __getitem__(self, i):
        if i == 0:
            return self.lo
        if i == 1:
            return self.hi
        else:
            raise IndexError("list index out of range")


if __name__ == "__main__":
    print("starting")
    i = Interval.from_tuple((1, 2))
    print(i)
    print(i[1])
    # print(Interval(0, 10).compare(Interval(2, 10)))
    # print(Interval(0.2, 10).compare(Interval(0.99, 10)))
    # s = Interval(0, 10); print(s)
    print(Interval(1676.162799346671, 1698).contains(Interval(1690, 1698)))
    print(Interval(1676.540526197659, 1711).contains(Interval(1690, 1699)))
