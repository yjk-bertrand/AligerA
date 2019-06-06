#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2018 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Scripts for running STEP1 of the AligerA pipeline tool.
STEP1 aims at aligning sequences using MAFFT.
The alignment step is not performed if the number of sequences N in the file exceeds the
[upper_sequence_limit] parameter.
If [MAFFT_upper_limit_addfragments] < N < [upper_sequence_limit], a slow but accurate
alignment is created out of the [MAFFT_upper_limit_addfragments] longest sequences and
then remaining sequences are added while trying to minimize the rearrangments to the 
initial alignment.
"""

# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
import sys
import shutil

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from aligera_scripts.utilities import run_subprocess


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def change_sequences_name(input_fasta, output_fasta):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    new_records = {}
    for record in records:
        if record.name[:3] == "_R_":
            new_name = record.name[3:]
            new_records[new_name] = SeqRecord(
                record.seq, name=new_name, id=record.id[3:], description=""
            )
        else:
            new_records[record.name] = record
    sorted_keys = sorted(list(new_records.keys()))
    sorted_records = [new_records[key] for key in sorted_keys]
    SeqIO.write(sorted_records, output_fasta, "fasta")


def run_MAFFT_small(fasta, cfg, **kargs):
    """
    Perform a regular MAFFT alignment when the number of sequences is below
    the [upper_sequence_limit] value.
    """
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    basename = fasta.split(in_suffix)[0]
    cmd = "mafft {mafft_settings} {input_fasta} > {output_fasta}".format(
        mafft_settings=cfg["MAFFT_parameters_small"],
        input_fasta=fasta,
        output_fasta=basename + "_temp_aligned.fasta",
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        template = "An exception of type {0} occurred when trying to run command {2}. \
    Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)
    change_sequences_name(basename + "_temp_aligned.fasta", basename + out_suffix)

    os.remove(basename + "_temp_aligned.fasta")
    new_fasta = basename + out_suffix
    try:
        shutil.move(new_fasta, out_folder)
    except Exception as ex:
        template = "An exception of type {0} occurred when trying to move fasta {2} \
        to directory {3}. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args, new_fasta, out_folder)
        print(message)
        sys.exit(message)
    s_i = "STEP1 Done with fasta: {}".format(fasta)
    return (s_i, "None")


def fragment_large_fasta(fasta, cfg, basename):
    """
    Fragment the initial file into two files:
        one with suffix '_core.fasta' that contains the longest sequences.
        one with suffix '_additional.fasta' that holds the remaining sequences.
    """
    max_nbr_seqs = cfg["MAFFT_upper_limit_addfragments"]
    records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    keys = records.keys()
    sorted_keys = sorted(keys, key=lambda x: len(records[x].seq), reverse=True)
    additional_set = sorted_keys[max_nbr_seqs:]
    core_set = sorted_keys[:max_nbr_seqs]
    SeqIO.write([records[x] for x in core_set], basename + "_core.fasta", "fasta")
    SeqIO.write(
        [records[x] for x in additional_set], basename + "_additional.fasta", "fasta"
    )
    return (basename + "_core.fasta", basename + "_additional.fasta")


def run_MAFFT_large(fasta, cfg, **kargs):
    """
    When [MAFFT_upper_limit_addfragments] < #sequences < [upper_sequence_limit]
        the sequence file is fragmented and only the [MAFFT_upper_limit_addfragments]
        are aligned with regular settings, and only then are the remaining sequences
        added.
    """
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    basename = fasta.split(in_suffix)[0]
    template = "An exception of type {0} occurred when trying to run command {2}. \
    Arguments:\n{1!r}"
    core_fasta, additional_fasta = fragment_large_fasta(fasta, cfg, basename)
    cmd = "mafft {mafft_settings} {input_fasta} > {output_fasta}".format(
        mafft_settings=cfg["MAFFT_parameters_small"],
        input_fasta=core_fasta,
        output_fasta=basename + "_core_aligned.fasta",
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)
    cmd = "mafft {mafft_settings} {additional_fasta} {input_fasta} >\
                {output_fasta}".format(
        mafft_settings=cfg["MAFFT_parameters_large"],
        input_fasta=basename + "_core_aligned.fasta",
        output_fasta=basename + "_temp_aligned.fasta",
        additional_fasta=additional_fasta,
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)
    change_sequences_name(basename + "_temp_aligned.fasta", basename + out_suffix)

    file_to_move = basename + out_suffix
    try:
        shutil.move(file_to_move, out_folder)
    except:
        print(
            "[Error] unable to move file {0} to folder {1}.".format(
                file_to_move, out_folder
            )
        )
    for x in [
        core_fasta,
        additional_fasta,
        basename + "_core_aligned.fasta",
        basename + "_temp_aligned.fasta",
    ]:
        try:
            os.remove(x)
        except:
            pass
    s_i = "STEP1 Done with fasta: {}".format(fasta)
    return (s_i, "None")
