# -*- coding: utf-8 -*-
"""
Created on Fri Apr 29 19:00:24 2022

@author: tol28
"""


import numpy as np
from collections import namedtuple
from functools import reduce


class Heap(object):
    """Min-heap.
    Objects that have been inserted using push can be retrieved in sorted
    order by repeated use of pop or pop_safe.
    """

    __slots__ = ["_root", "_nelems"]

    def __init__(self, items=()):
        self._root = None
        self._nelems = 0
        # print("Initializing Heap...")
        for x in items:
            self.push(x)

    def __iadd__(self, other):
        """Merge other into self, destroying other in the process."""

        self._root = _meld(self._root, other._root)
        self._nelems += other._nelems

        other._root = None
        other._nelems = 0

        return self

    def __len__(self):
        return self._nelems

    def _pop(self):
        r = self._root.key
        self._root = _pair(self._root.sub)
        self._nelems -= 1
        return r

    def pop(self):
        """Remove the smallest element from the heap and return it.
        Raises IndexError when the heap is empty.
        """
        try:
            return self._pop()
        except AttributeError:
            raise IndexError("pop from an empty Heap")

    def pop_safe(self):
        """Like pop, but returns None when the heap is empty."""
        return self._root and self._pop()

    def push(self, x):
        """Push element x onto the heap."""
        self._root = _meld(self._root, _Node(x, []))
        self._nelems += 1

    @property
    def top(self):
        """The smallest element of the heap."""
        try:
            return self._root.key
        except AttributeError:
            raise IndexError("min of an empty Heap")


_Node = namedtuple("_Node", "key sub")


def _meld(l, r):
    """Meld (merge) two pairing heaps, destructively."""
    # We deviate from the usual (persistent) treatment of pairing heaps by
    # using list's destructive, amortized O(1) append rather than a "cons".
    if l is None:
        return r
    elif r is None:
        return l
    elif l.key[-1] < r.key[-1]:
        l.sub.append(r)
        return l
    else:
        r.sub.append(l)
        return r


def _pair(heaps):
    """Pair up (recursively meld) a list of heaps."""
    return reduce(_meld, heaps, None)


def heap_dijkstra_(G,c,source,flag_return_pred=False):
    H = Heap()
    dist = {}
    for j in G.nodes():
        dist[j] = np.float("inf")
    dist[source] = 0
    pred = {source : 0}
    H.push((source,0))
    tabu = []
    while H.__len__()>0:
        i = H.pop()
        if i[0] in tabu: continue
        tabu.append(i[0])
        for e in G.edges(i[0]):
            j = e[0] if e[0]!=i[0] else e[1]
            value = dist[i[0]]
            if G.is_directed():
                value += c[e]
            else:
                try:
                    value += c[e]
                except KeyError:
                    value += c[(e[1],e[0])]
            if dist[j]>value:
                pred[j] = i[0]
                H.push((j,value))    
                dist[j] = value
    if not flag_return_pred:
        return dist
    else:
        return dist,pred
