#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import os

# from memory_profiler import profile
import shelve
import time
from multiprocessing import Manager
import concurrent.futures


# =======================================================================================
#               CLASS
# =======================================================================================


class UnionFind:
    """Union-find data structure with path compression and Union-by-rank.
    getComponents return a list of sets, each set instance represents a 
    component.  
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}
        self.components = {}

    def __getitem__(self, object):
        """
        Find and return the name of the set containing the object.
        """
        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            self.components[object] = [object]
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]
        for ancestor in path:
            self.parents[ancestor] = root
            if ancestor != root:
                self.components[root].append(ancestor)
                try:
                    del self.components[ancestor]
                except:
                    continue
        return root

    def union(self, objects):
        """
        Find the sets containing the objects and merge them all.
        """
        roots = [self[x] for x in objects]  # uses__getitem__
        heaviest = max([(self.weights[r], r) for r in roots])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest
                self.components[heaviest].extend(self.components[r])
                del self.components[r]

    def getComponents(self):
        return [set(x) for x in self.components.values()]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


def grouper(n, iterable, padvalue=None):
    return [iterable[i::n] for i in range(n)]


def _sequential(k, dict_elts, shl):
    keys = dict_elts[k]
    with UnionFind() as UV:
        for i in keys:
            for item in shl[i]:
                UV.union(item)
    return k, UV.getComponents()


def group_tasks(queue, result_dict):
    try:
        q = queue.get(True, 0.05)
        result_dict[q] = None
    except queue.Empty:
        print("Nothing to be done, queue is empty")


def process_future(script_to_run, list_elts, result_dict, queue, cpu_number, shl):
    """
    compute components graph using multiprocessing and ProcessPoolExecutor
    list_elts is a list of partitions of shelve keys
    result_dict is the dictionary linked to the Manager
    shl is an open shelve with species as key and graph edges as values
    """
    dict_elts = dict(enumerate(list_elts))
    keys = list(dict_elts.keys())
    for k in keys:
        queue.put(k)

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=cpu_number
    ) as group_processes:
        for i in range(queue.qsize()):
            group_processes.submit(group_tasks, queue, result_dict)

    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_number) as processes:
        future_tasks = {
            processes.submit(script_to_run, k, dict_elts, shl): k
            for k in range(cpu_number)
        }

        for future in concurrent.futures.as_completed(future_tasks):
            result_dict.append(future.result()[1])


def _parallel(shl, n_partitions):
    sendbuf = grouper(n_partitions, list(shl.keys()))
    manager = Manager()
    queue = manager.Queue()
    result_dict = manager.list()
    process_future(_sequential, sendbuf, result_dict, queue, n_partitions, shl)
    return list(result_dict)


if __name__ == "__main__":
    print("starting")
    # start_time = time.time()
    os.chdir(r"E:\\python\\aligera\\aligera_v_1_0_b\\temp")
    shl = shelve.open("common_hits_dict")
    start_time = time.time()
    # components = _sequential_non_parallel(shl)
    # print(list(components))
    components = _parallel(shl, 4)
    print("Number of components: ", len(components))
    print([len(c) for c in components])
    print(len(merge(components)))
    print(" components extraction in {} seconds".format(time.time() - start_time))
    shl.close()
