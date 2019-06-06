#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
#   Class Node which contains a data structure use in the interval search tree.
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

try:
    from aligera_scripts.interval import Interval
except:
    from interval import Interval


# =======================================================================================
#               CLASS
# =======================================================================================


class Node:
    """
    Data structure use in the interval search tree. It relies on the node data 
    structure found at:
    https://algs4.cs.princeton.edu/93intersection/IntervalST.java.html
    This data structure is used in the STEP2 of the AligerA pipeline tool to trim
    the flanking regions of an alignment.
    """

    def __init__(self, interval, value, left_node=None, right_node=None):
        self.interval = interval
        self.value = [value]
        self.left = left_node
        self.right = right_node
        self.N = 1  # size of subtree
        self.max = interval.getMax  # max endpoint in subtree rooted at this note


    @classmethod
    def from_tuple(cls, tuple_lo_hi, value):
        if len(tuple_lo_hi) != 2:
            raise ValueError(
                "Node: tuple can contain only two values {} ".format(tuple_lo_hi)
            )
        node = cls(Interval(tuple_lo_hi[0], tuple_lo_hi[1]), value)
        return node


    def __str__(self):
        return "Node with interval: {0} and value {1}".format(self.interval, self.value)


if __name__ == "__main__":
    print("starting")
    n1 = Node.from_tuple((1, 2), 2)
    # n1.value.extend([3])
    print(n1)
