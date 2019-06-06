# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
#   Merge lists of components.

__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import itertools
import sys


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def mergeComponents(L):
    """
    Merge lists of components.
    L is a list that contains lists of components
    Return a bag that for each components yields the components that intersects.
    """
    if not L:
        s = "[Error] empty list no components to merge: exiting"
        sys.exit(s)
    if len(L) == 1:
        return [[] for x in range(len(L))]
    L0 = L[0]
    L1 = L[1]
    M = [[] for x in range(len(L0) + len(L1))]
    x = 0

    for i in L0:
        y = len(L0)
        for j in L1:
            if hasMerged(M, x, y):
                M[x].append(y)
                M[y].append(x)
            elif not i.isdisjoint(j):
                M[x].append(y)
                M[y].append(x)
            y += 1
        x += 1
    return M


def hasMerged(M, s, t):
    """
    Helper function for merge_components() 
    Instead of checking for component intersection directly it analyses 
    the result bag to find transitive intersections
    s = source in L0, t = target in L1
    """
    node_s = set(M[s])
    for n in M[t]:
        if not node_s.isdisjoint(set(M[n])):
            return True
    return False


def getMergedComponents(C, M):
    """
    translate the merging bag M into a new components list
    C is the original list of components and M is the bag of components
    to merge produced by mergeComponents()
    """
    results = []
    used_components = [0 for x in range(len(M))]
    L = list(itertools.chain.from_iterable(C))
    for idx, elt in enumerate(M):
        #  check if components has been consumed
        if used_components[idx]:
            continue
        if not elt:
            results.append(L[idx])
        else:
            components = [L[idx]]
            used_components[idx] = 1
            queue = elt
            while queue:
                i = queue.pop()
                if not used_components[i]:
                    used_components[i] = 1
                    components.append(L[i])
                    queue.extend(M[i])
            results.append(set.union(*components))
    return results


def merge(L_components):
    """
    Merge lists of components.
    """
    while len(L_components) > 1:
        L0 = L_components.pop()
        L1 = L_components.pop()
        M = mergeComponents([L0, L1])
        L_components.append(getMergedComponents([L0, L1], M))
    return L_components[0]


if __name__ == "__main__":
    print("starting")
    components_nodes = [
        [set([1, 2, 3, 4]), set([6, 7, 8]), set([9, 10])],
        [set([4, 5]), set([6, 7, 8, 9]), set([11])],
        [set([1, 2, 3, 4]), set([6, 7, 8]), set([12])],
        [set([1, 2, 3, 4]), set([6, 7]), set([13, 8])],
    ]
    result_components = merge(components_nodes)
    print(result_components)

    print(
        getMergedComponents(
            [[{1, 2, 3, 4}, {6, 7}, {8, 13}], [{1, 2, 3, 4}, {8, 6, 7}, {12}]],
            [[3], [4], [4], [0], [1, 2], []],
        )
    )
