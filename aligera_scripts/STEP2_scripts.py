#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
    Functions used by the AligerA pipeline tool to trim an alignment out of sparsely
    overlapping flanking regions and highly gapped columns.
    Finding flanking limits uses an interval search tree algorithm and removing 
    mostly empty columns uses a sliding window.
"""
__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import shutil

try:
    from aligera_scripts.find_limits import FindLimits
except:
    from find_limits import FindLimits

from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def find_gapped_columns(align, cfg):
    """Find all columns that contain more gaps than the cfg setting using 
    a sliding windows"""
    max_gap_proportion = cfg["max_gap_proportion"]
    nbr_sequences = len(align)
    columns_to_remove = []
    len_align = align.get_alignment_length()
    for index in range(len_align):
        column = str(align[:, index]).replace("n", "N")
        gap_freq = (column.count("N") + column.count("-")) / nbr_sequences
        if gap_freq > max_gap_proportion:
            columns_to_remove.append(index)
    if columns_to_remove:
        idxs = [x for x in range(len_align) if x not in columns_to_remove]
        trimmed_records = []
        for rec in align:
            L_seq = list(rec.seq)
            new_seq = "".join([L_seq[i] for i in idxs])
            new_rec = SeqRecord(
                Seq(new_seq, IUPAC.IUPACAmbiguousDNA()),
                name=rec.name,
                id=rec.id,
                description="",
            )
            trimmed_records.append(new_rec)
        return trimmed_records
    return align


def trim(fasta, cfg, **kargs):
    """
    Main function for STEP2:
    interp. 
        ungap_proportion (float)    Minimum proportion of non gap positions to 
                                    be considered during limit search.
        min_aa_proportion (float)   Minimum proportion of conserved aa in the column
                                    for the column to be considered during limit search.
        max_gap_proportion (float)  Maximum proportion of gaps in the column. 
                                    Columns that have higher proportions of gap are 
                                    removed.
        window_length (int)         length of the window in nucleotides that is 
                                    translated in order to search for conserved 
                                    consecutive columns of nucleotides.                                     
    """
    window_length = (cfg["window_length"] // 3) * 3
    if window_length < 3:
        exception = "[Error] Problem window length < 3 \n\
        Exiting..."
        raise Exception(exception)
        sys.exit()
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    basename = fasta.split(in_suffix)[0]

    try:
        assert (
            0 <= cfg["ungap_proportion"] <= 1
            and 0 <= cfg["min_aa_proportion"] <= 1
            and 0 <= cfg["max_gap_proportion"] <= 1
        )
    except:
        s = "Wrong parameters: 'ungap_proportion', 'min_aa_proportion', \
        'max_gap_proportion' must be in [0, 1]"
        raise ValueError(s)
        sys.exit(s)
    try:
        align = AlignIO.read(fasta, "fasta")
    except:
        exception = "[Error] Problem opening fasta alignment {}\n\
        Are all sequences the same length?".format(
            fasta
        )
        raise Exception(exception)
        return exception

    n_seq = len(align)
    if n_seq < 3:
        s_i = "STEP2 Done with fasta: {}".format(fasta)
        s_d = "No transformation performed on fasta {0}\n\
        which only contains {1} sequences".format(
            fasta, n_seq
        )
        return (s_i, s_d)

    try:
        FL = FindLimits(align, **cfg)
    except Exception as ex:
        template = "An exception of type {0} occurred when trying to\
 instantiate FindLimit object. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        sys.exit(0)

    FL.findFragments()
    FL.findLimits()
    align_length = align.get_alignment_length()
    start = FL.start
    end = FL.end
    if start > end:
        exception = "[Error] Problem with finding start/end \
        for fasta {0}; Start: {1}, End: {2}".format(
            fasta, start, end
        )
        raise Exception(exception)
        start = 0
        end = align_length

    new_align = align[:, start:end]
    trimmed_align = find_gapped_columns(new_align, cfg)
    SeqIO.write(trimmed_align, basename + out_suffix, "fasta")
    try:
        shutil.move(basename + out_suffix, out_folder)
    except:
        exception = "[Error] Unable to move file {0} to folder {1}".format(
            basename + out_suffix, out_folder
        )
        raise Exception(exception)

    s_i = "STEP2 Done with fasta: {}".format(fasta)
    s_d = "start: {0} end: {1}".format(start, end)
    return (s_i, s_d)

