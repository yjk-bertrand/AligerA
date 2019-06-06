#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
#  Class IntervalST that create a interval search tree used to find the limits of
#  an alignment

__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================
try:
    from aligera_scripts.node import Node
except:
    from node import Node

try:
    from aligera_scripts.interval import Interval
except:
    from interval import Interval


# =======================================================================================
#               CLASS
# =======================================================================================


class IntervalST:
    """
    Create interval search tree using a modified version of the IntervalST data object:
        https://algs4.cs.princeton.edu/93intersection/IntervalST.java.html
        This data structure is used in the STEP2 of the AligerA pipeline tool to trim
        the flanking regions of an alignment.
    """

    def __init__(self, root=None):
        self.root = root


    def put(self, interval, value):
        """
        New data object in inserted at the root of the tree (no random insertion). 
        Value is a list.
        If interval is alredy present in the tree, the new value is appended to the
        value list of the interval.
        """
        self.root = self._rootInsert(self.root, interval, value)


    def _rootInsert(self, node, interval, value):
        if node is None:
            return Node(interval, value)
        cmp = interval.compare(node.interval)
        if cmp < 0: #  adding left
            node.left = self._rootInsert(node.left, interval, value)
            node = self._rotR(node)
        if cmp > 0: #  adding right
            node.right = self._rootInsert(node.right, interval, value)
            node = self._rotL(node)
        else:
            if value not in node.value:
                node.value.append(value)
        return node


    def get(self, interval):
        """return value associated with the given key
            if no such value, return null"""
        return self._get(self.root, interval)


    def _get(self, node, interval):
        if node is None:
            return None
        cmp = interval.compare(node.interval)
        if cmp < 0: #  searching left
            node.left = self._get(node.left, interval)
        elif cmp > 0:  #  searching right
            node.right = self._get(node.right, interval)
        else:
            return node.value


    def _rotR(self, node):
        """
        Rotate tree right
        """
        x = node.left
        node.left = x.right
        x.right = node
        self._fix(node)
        self._fix(x)
        return x


    def _rotL(self, node):
        """
        Rotate tree left
        """
        x = node.right
        node.right = x.left
        x.left = node
        self._fix(node)
        self._fix(x)
        return x


    def _fix(self, node):
        """
        fix auxilliar information (subtree count and max fields)
        """
        if node is None:
            return
        node.N = 1 + self._size(node.left) + self._size(node.right)
        node.max = max(
            node.interval.hi, self._nodeMax(node.left), self._nodeMax(node.right)
        )


    def _nodeMax(self, node):
        if node == None:
            return 0
        else:
            return node.max


    def size(self):
        return self._size(self.root)


    def searchContains(self, interval):
        """
        Main algorithm for searching the interval tree. In comparison to the classical
        algorithm it has been modified to work with included interval instead of mere
        overlapping intervals.
        """
        L = []
        self._searchContains(self.root, interval, L)
        return L


    def _searchContains(self, node, interval, L):
        found1 = False
        found2 = False
        found3 = False
        if node is None:
            return False
        if node.interval.contains(interval):
            L.append(node)
            found1 = True
        if (
            node.left is not None
            and node.left.max >= interval.hi
            and node.left.max >= interval.lo
        ):
            found2 = self._searchContains(node.left, interval, L)
        if (node.right is not None) and (
            found2
            or (node.left is None)
            or node.left.max <= interval.lo
            or (node.right.max >= interval.hi and node.right.max >= interval.lo)
        ):
            found3 = self._searchContains(node.right, interval, L)

        return found1 or found2 or found3


    def searchIntersect(self, interval, value):
        """
        Algorithm for searching for intersecting intervals that contain a different
        value to the one in the arguments. Return a boolean.
        """
        return self._searchIntersect(self.root, interval, value)


    def _searchIntersect(self, node, interval, value):
        while node != None:
            if (
                interval.intersects(node.interval)
                and len(set(node.value).difference(set([value]))) > 1
            ):
                return True
            elif node.left == None:
                node = node.right
            elif node.left.max < interval.lo:
                node = node.right
            else:
                node = node.left
        return False


    def height(self):
        """
        Height of the Search Tree
        """
        return self._height(self.root)


    def _height(self, node):
        if node is None:
            return 0
        return 1 + max(self._height(node.left), self._height(node.right))


    def _size(self, node):
        if node == None:
            return 0
        else:
            return node.N


if __name__ == "__main__":
    print("starting")
    ST = IntervalST()
    ST.put(Interval(5, 10), 1)
    ST.put(Interval(6, 11), 2)
    ST.put(Interval(6, 11), 3)
    ST.put(Interval(20, 25), 7)
    print(ST.searchIntersect(Interval(20, 25), 7))

    # print("result value is ", ST.get(Interval(1, 6)))
    # print("search:", [x.value for x in ST.search(Interval(20, 21))])
#    import pickle
#    ST = pickle.load(open("ST.dat", 'rb'))
#    print("root is ", ST.root)
#    node = ST.root.right.right.left.right
#    print(node.left)
#    print(ST.root.max)
#    print("height is:", ST.height())

#    inter =  pickle.load(open("intervals.dat", 'rb'))
#    for i in inter:#[12:14]:
#        print("search interval", i)
#        l = [print(x) for x in ST.search(i)]
#        print("interval length", len(l))
#        print("\n")
#        if len(l) > 2: break

#    print("tree size:", ST.size())
