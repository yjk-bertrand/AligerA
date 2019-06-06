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
import sys
import math
from collections import Counter

from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC, Gapped
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment as MulAlign

try:
    from aligera_scripts.interval import Interval
except:
    from interval import Interval
try:
    from aligera_scripts.interval_st import IntervalST
except:
    from interval_st import IntervalST


# =======================================================================================
#               CLASS
# =======================================================================================


class FindLimits:
    """
    Algorithn used in AligerA in the STEP2 of the pipeline. It trims the flanking 
    regions based on the first occurance on each side of the alignment of a conserved
    column of amino acids. The algorithm relies on an interval search tree data
    structure to efficiently determine the overlapping regions of the alignment.
    """

    def __init__(self, alignment, **kargs):
        """Constructor:
            alignment: a Bio alignment object
            window_length: number of nucleotides that are translated during
                window slide
            ungap_ratio: ratio of (non gap elt in column)/(column length)
                column is used in translation when ratio is higher than set value
            min_aa_ratio: ratio of (largest count of aa)/(total number of aa)
                column is used in limit search when ratio is higher
        """
        self.align = alignment
        self.length = alignment.get_alignment_length()
        #  Make sure the windows length is a multiple of 3
        self.window_length = (kargs["window_length"] // 3) * 3
        #  Compute the number of sequences to use out of a proportion value
        self.min_number_seq = math.ceil(kargs["ungap_proportion"] * len(alignment))
        if self.min_number_seq < 1:
            raise ValueError(
                '[Error] Too few sequences or "ungap_proportion" set too low'
            )
        self.min_aa_ratio = kargs["min_aa_proportion"]
        if self.window_length < 3:
            raise ValueError("[Error] Window length too short (<3)")
        #  Sequence fragments longer than window length
        self.fragments = []
        #  Sequence fragments longer than window length that intersect at least one
        #  more fragment
        self.fragmentsIntersect = []
        #  Three translation frames
        self.t_align0 = self.translate(self.align, 0)
        self.t_align1 = self.translate(self.align, 1)
        self.t_align2 = self.translate(self.align, 2)
        self.start = 0
        self.end = self.length


    def findFragments(self):
        """
        For all sequences find fragments that are longer than window_length parameter.
        """
        align_length = self.align.get_alignment_length()
        idx = -1
        for rec in list(self.align):
            idx += 1
            seq = rec.seq
            f = []
            i = 0
            start = 0
            while i < align_length:
                while i < align_length and seq[i] in ["-", "N", "n"]:
                    i += 1
                start = i
                while i < align_length and seq[i] not in ["-", "N", "n"]:
                    i += 1
                end = i
                if (end - start) >= self.window_length:
                    f.append(((start, end), idx))
            self.fragments.extend(f)


    def findLimits(self):
        """
        Use the fragments limits self.fragments to build an interval search tree and 
        search through the fragments to find overlapping conserved regions.
        """
        try:
            assert len(self.fragments) > 0
        except:
            exception = "Fragment list empty, no limits to be extracted"
            print(exception)
            return

        ST = IntervalST()

        for (limits, idx) in self.fragments:
            ST.put(Interval.from_tuple(limits), idx)

        for (limits, idx) in self.fragments:
            if ST.searchIntersect(Interval.from_tuple(limits), idx):
                self.fragmentsIntersect.append((limits, idx))

        try:
            assert len(self.fragmentsIntersect) > 0
        except:
            exception = "No fragment intersect any other, no limits to be extracted"
            print(exception)
            return

        try:
            self.getStart(ST)
        except Exception as ex:
            template = "An exception of type {0} occurred when trying to \
        find lower limit. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)
            sys.exit(0)

        try:
            self.getEnd(ST)
        except Exception as ex:
            template = "An exception of type {0} occurred when trying to \
        find upper limit. Arguments:\n{1!r}"
            message = template.format(type(ex).__name__, ex.args)
            print(message)
            sys.exit(0)


    def getStart(self, interval_ST):
        """
        Infer start
        """
        try:
            assert len(self.fragmentsIntersect) > 0 and interval_ST.size() > 0
        except:
            exception = "Fragment list empty or interval search tree of null \
        size no lower limit to infer"
            print(exception)
            return

        visited = [0 for x in range(self.length)]
        for fragment in sorted(self.fragmentsIntersect):
            #  Work on each fragment
            (b, e), idx = fragment
            begin = math.floor(b)
            l = ((e - begin) // 3) * 3
            end = begin + l
            triplets = [
                (i, i + self.window_length)
                for i in range(begin, end - self.window_length)
            ]
            intervals = [Interval.from_tuple(t) for t in triplets]
            for i in intervals:
                idx = i[0]
                if visited[idx]:
                    continue
                visited[idx] = 1
                intersect = []
                for x in interval_ST.searchContains(i):
                    intersect.extend(x.value)
                if len(intersect) >= self.min_number_seq:
                    results = self.testLimit(intersect, begin)
                    if all(results):
                        self.start = i[0]
                        return


    def getEnd(self, interval_ST):
        """
        Infer end
        """
        try:
            assert len(self.fragmentsIntersect) > 0 and interval_ST.size() > 0
        except:
            exception = "Fragment list empty or interval search tree of null \
            size no upper limit to infer"
            print(exception)
            return

        visited = [0 for x in range(self.length)]
        for fragment in sorted(
            self.fragmentsIntersect, key=lambda x: x[0][1], reverse=True
        ):
            #  Work on each fragment
            (b, end), idx = fragment
            l = ((end - math.floor(b)) // 3) * 3
            begin = end - l 
            triplets = [
                (i, i + self.window_length)
                for i in reversed(range(begin, end - self.window_length + 1))
            ]
            intervals = [Interval.from_tuple(t) for t in triplets]
            for i in intervals:
                intersect = []
                idx = i[0]
                if visited[idx]:
                    continue
                visited[idx] = 1
                for x in interval_ST.searchContains(i):
                    intersect.extend(x.value)
                if len(intersect) >= self.min_number_seq:
                    results = self.testLimit(intersect, begin)
                    if all(results):
                        self.end = i[1]
                        return


    def testLimit(self, list_seqs, start):
        """
        Extract the aa sequences in the window.
        list_seqs is the list of sequence id in the alignment (not the id associated
            with the Bio.Seq object).
        start is the index of the start of the window.
        """
        frame = start % 3
        aa_window_length = int(self.window_length / 3)
        begin = int((start - frame) / 3)
        end = int(begin + aa_window_length)

        if frame == 0:
            t_align = self.t_align0
        elif frame == 1:
            t_align = self.t_align1
        else:
            t_align = self.t_align2

        sub_align = MulAlign([], Gapped(IUPAC.ExtendedIUPACProtein(), "N"))
        for idx in list_seqs:
            sub_align.append(t_align[idx][begin:end])

        result = []

        for c in range(aa_window_length):
            c = Counter(sub_align[:, c])
            #  count the most common aa
            nbr_most_common = c.most_common(1)[0][1]
            if nbr_most_common / len(list_seqs) >= self.min_aa_ratio:
                result.append(True)
            else:
                result.append(False)
        return result


    def translate(self, align, offset):
        """
        Translate the alignment according to the selected frame which is set 
            according to 'offset' value
        """
        end = ((align.get_alignment_length() - offset) // 3) * 3 + offset
        t_align = MulAlign([], Gapped(IUPAC.ExtendedIUPACProtein(), "N"))
        for rec in align:
            seq = str(rec.seq).upper().replace("-", "N").replace("n", "N")
            new_seq = Seq(seq, IUPAC.IUPACAmbiguousDNA())[offset:end].translate()
            new_rec = SeqRecord(new_seq, name=rec.name, id=rec.id, description="")
            t_align.append(new_rec)
        return t_align



