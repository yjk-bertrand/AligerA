#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Clean all sequence regions that have low similarity to the other sequences in 
the alignment.
The alignments are filtered for sequences that are too short or contain too many
ambiguous positions. Ambiguous positions hold any of the IUPAC degeneracy symbols.
Taxa are infered from sequences names and alignments that contain a number of taxa
above the selected threshold are saved with the suffix '_complete.fasta' otherwise
they get the suffix '_incomplete.fasta'.
"""
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import shutil
import itertools
import math
from collections import Counter

try:
    from aligera_scripts.utilities import get_alleles
except:
    from utilities import get_alleles

from Bio import SeqIO


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def calculated_hamming_distance_kmer(kmer_1, kmer_2):
    """
    Calculate Hamming distance for pair of kmers
    """
    return sum(a == b for a, b in zip(kmer_1, kmer_2) if a != "-" and b != "-")


def find_singletons_per_position(taxa, aln_length, kmer, sensitivity, alleles):
    """
    interp. taxa (dict)         taxa number as key, sequences as values
            allele (list(set))  Allele pairs are in sets
            
    A dictionary of dictionaries is generated
    D.items() >>>[(2, {1: '--T'}),  (6, {0: '--A', 1: 'ATT'})
    with keys as position and within subdictionaries, taxa as keys, sequences as values.   
    """
    aln_dict = {}
    for name in list(taxa.keys()):
        #  For each sequence create a list of (index, kmer) starting at each position.
        #  The aln_dict is populated with positions as keys and dict as values
        #   where seq names are keys and kmers are values.
        s_kmers = list(
            enumerate([taxa[name][idx : idx + kmer] for idx in range(aln_length)])
        )
        for item in s_kmers:
            if item[1] != "-" * kmer:
                if item[0] in list(aln_dict.keys()):
                    aln_dict[item[0]][name] = item[1]
                else:
                    aln_dict[item[0]] = {name: item[1]}

    aln_dict_no_duplicates = {}
    #  Yields a dictionary (k, v) with position as k and set of kmers at that position:
    #  {2: ['--T'], 3: ['-TA'], 4: ['TTC', '-AC'], 7: ['ACA', 'TCA']}
    for k, v in aln_dict.items():
        aln_dict_no_duplicates[k] = list(set(v.values()))

    counter_dict = {}
    #  Dictionary that hold for each position a kmer counter
    for k, v in aln_dict.items():
        counter_dict[k] = Counter(v.values())

    #  Dictionary that for each position (key) holds kmer as key and
    #  as values, a list of sequences that possess that kmer:
    #   {2: {'--T': [1, 2], '-CG':[0, 3]}, 3: {'-TA': [0,1,2,3]}}
    #  If two positions have the same kmer the first is overwritten.
    aln_dict_reversed = {x: {} for x in aln_dict.keys()}
    for idx in aln_dict.keys():
        D = aln_dict_reversed[idx]
        for k, v in aln_dict[idx].items():
            if v not in D.keys():
                D[v] = [k]
            else:
                D[v].append(k)
    # list of positions made only of ambiguous bases
    non_singletons_per_position_list = []

    #  For all non empty position compute hamming distances (HD) between pairs of Kmer.
    #  If HD > sensitivity, accept pairs. Sets of intersecting pairs are created
    for idx in range(aln_length):
        if idx not in list(aln_dict_no_duplicates.keys()):
            non_singletons_per_position_list.append([])
        else:
            kmers_list = aln_dict_no_duplicates[idx]
            L_idxs_intermediate = set()
            for k_pair in list(itertools.combinations(kmers_list, 2)):
                HD = calculated_hamming_distance_kmer(k_pair[0], k_pair[1])
                if HD > sensitivity:
                    edge = set([k_pair[0], k_pair[1]])
                    if not L_idxs_intermediate:
                        L_idxs_intermediate = edge
                    else:
                        L_idxs_intermediate = edge.union(L_idxs_intermediate)
            if L_idxs_intermediate:
                non_singletons_per_position_list.append(L_idxs_intermediate)
            else:
                non_singletons_per_position_list.append([])

    singletons_per_position_list = []
    #  For each position (idx), we already have the non singletons kmer.
    #  The potential singleton kmer are the kmers that are not present in at least
    #  two sequences within the Hamming distance limit.
    #  For each potential singleton we check that it really occurs
    #  only once with the counter, then we extract the corresponding kmer.
    #  If it occurs twice, add the pair of sequence to the list of sequence to correct
    #  if the pair is an allelic pair.
    for idx, st in enumerate(non_singletons_per_position_list):
        if idx not in list(aln_dict_no_duplicates.keys()):
            singletons_per_position_list.append([])
        else:
            not_present = list(set(aln_dict_no_duplicates[idx]).difference(st))
            if not not_present:
                singletons_per_position_list.append([])
            else:
                true_singletons = []
                for kmer in not_present:
                    if counter_dict[idx][kmer] == 1 or (
                        counter_dict[idx][kmer] == 2
                        and set(aln_dict_reversed[idx][kmer]) in alleles
                    ):
                        true_singletons.extend(aln_dict_reversed[idx][kmer])
                singletons_per_position_list.append(true_singletons)
    return singletons_per_position_list


def find_start_end(seq, kmer):
    """
    Retrieve the start and the end of the sequence,
    i.e. the first position that is followed by kmer nucleotides
    and the last position that end with kmer nucleotides
    """
    for idx in range(len(seq)):
        start = 0
        if seq[idx] != "-" and "-" not in seq[idx : idx + kmer]:
            start = idx
            break
    for idx in range(len(seq)):
        end = len(seq)
        if seq[-idx] != "-" and "-" not in reversed(seq[-idx : -idx - kmer]):
            end = len(seq) - idx
            break
    return (start, end)


def make_split_dict(singletons_per_position, inverted_taxa_dict):
    """
    Create a dictionary that record the splits and their positions
    """
    split_dict = {}
    for idx in range(len(singletons_per_position)):
        for taxon in [
            inverted_taxa_dict[x] for x in list(singletons_per_position[idx])
        ]:
            if taxon not in list(split_dict.keys()):
                split_dict[taxon] = [idx]
            else:
                split_dict[taxon].append(idx)
    return split_dict


def write_fasta(fasta_name, sequence_dict):
    with open(fasta_name, "w") as f:
        for taxon in sorted(list(sequence_dict.keys())):
            f.write(">{}\n".format(taxon))
            f.write("{}\n".format(sequence_dict[taxon]))


def check_taxa_name_format(seq_name, fasta):
    """
    Check whether the sequence name
    has been properly formated, with '|' separating
    the taxon name from the sequence id
    """
    if "|" not in seq_name:
        exception = "[Error] Format problem with sequence {0} in fasta {1}\n\
        the sequence name should contain a '|' separation between \n\
        the taxon name and the sequence id".format(
            seq_name, fasta
        )
        raise Exception(exception)


def get_seq(record):
    """
    Return the sequence as a string replacing 'N'/'n' with '-'
    """
    seq = str(record.seq).replace("N", "-").replace("n", "-")
    return seq


def find_unspliced_segments(fasta, cfg, **kargs):
    """
    Main function for STEP3:
    interp.
        kmer (int)                      The number of consecutive bases to be examine
                                        for computing Hamming distance.                                             
        max_Hamming_distance (int)      Max value of dissimilarity between two kmers 
                                        for the Hamming distance (HD) to be 1.
                                        max_Hamming_distance < kmer.
                                        E.g. kmer = 10 max_Hamming_distance = 6, 
                                        if similarity is <= 6 , HD = 0, 
                                        if similarity > 6, HD = 1.
        min_length_ambiguous_sequence (int) Minimum length of the window for searching
                                        for faulty sequence region.
        consecutive_ambiguity (int)     Maximum number of 'N' allowed in the window to
                                        be tagged as ambiguous. 
                                        consecutive_ambiguity < min_length_ambiguous_sequence
        difference_proportion (int)     [0,1] Minimum proportion of positions in the 
                                        window that have HD = 0, in order to indentify 
                                        a faulty sequence region. 
        min_taxa_in_alignment (int)     > 0 Alignments with more taxa are saved with
                                        the suffix '_complete.fasta' otherwise 
                                        they get the suffix '_incomplete.fasta'.
        outgroups (List of str)         Taxa in the outgroup are not processed.
        min_alignment_length (int)      After removing faulty regions, a sequence is 
                                        kept in the alignment if its length is above 
                                        this value.
    """
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    kmer = cfg["kmer"]
    outgroups = cfg["outgroups"]
    min_aln_length = cfg["min_alignment_length"]
    sensitivity = cfg["max_Hamming_distance"]
    ambiguities = cfg["consecutive_ambiguity"]
    max_ambiguity = cfg["max_ambiguity_proportion"]
    min_taxa = cfg["min_taxa_in_alignment"]
    min_length_fragment = cfg["min_length_ambiguous_sequence"]
    differences = math.floor(min_length_fragment * cfg["difference_proportion"])
    basename = fasta.split(in_suffix)[0]
    error_template = "An exception of type {0} occurred when trying to \
 {2}. Arguments:\n{1!r}"
    try:
        assert (min_taxa) > 0
    except:
        s = "Wrong parameters: 'min_taxa_in_alignment' must be > 1"
        raise ValueError(s)
        sys.exit(s)

    try:
        assert min_aln_length > min_length_fragment
    except:
        s = "Wrong parameters: 'min_alignment_length' must be > \
        min_length_ambiguous_sequence"
        raise ValueError(s)
        sys.exit(s)

    try:
        assert cfg["difference_proportion"] <= 1
    except:
        s = "Wrong parameters: 'difference_proportion' must be <= 1"
        raise ValueError(s)
        sys.exit(s)

    try:
        assert sensitivity <= kmer
    except:
        s = "Wrong parameters: 'max_Hamming_distance' must be <= 'kmer'"
        raise ValueError(s)
        sys.exit(s)

    try:
        assert ambiguities < min_length_fragment
    except:
        s = "Wrong parameters: min_length_ambiguous_sequence must be longer \
        than 'consecutive_ambiguity'"
        raise ValueError(s)
        sys.exit(s)

    try:
        records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    except Exception as ex:
        s = "open alignment: {}".format(fasta)
        message = error_template.format(type(ex).__name__, ex.args, s)
        print(message)
        sys.exit(0)
    n_seq = len(records)
    if n_seq < 3:
        s_i = "STEP3 Done with fasta: {}".format(fasta)
        s_d = "No transformation performed on fasta {0}\n\
        which only contains {1} sequences".format(
            fasta, n_seq
        )
        return (s_i, s_d)

    #  1- Remove empty sequences
    remove_key = []
    for rec in records.keys():
        check_taxa_name_format(records[rec].name, fasta)
        if set(list(get_seq(records[rec]))) == set(["-"]):
            remove_key.append(rec)

    for rec in remove_key:
        records.pop(rec, None)

    #  2- Filter out sequences from the outgroup that might be too divergent
    #     and therefore can be mistaken for misalignments
    #     taxa_dict holds the taxa names in order to free memory
    
    
    taxa_names = list(
        enumerate(
            [
                records[rec].name
                for rec in records.keys()
                if rec.split("|")[0] not in outgroups
            ]
        )
    )

    taxa_dict = {x[1]: x[0] for x in taxa_names}  #  dict name:number
    inverted_taxa_dict = {x[0]: x[1] for x in taxa_names}  #  dict number:name
    taxa = {
        taxa_dict[records[rec].name]: get_seq(records[rec])
        for rec in records.keys()
        if rec.split("|")[0] not in outgroups
    }

    alleles = [
        set((taxa_dict[x], taxa_dict[y]))
        for (x, y) in get_alleles(records.keys(), outgroups)
    ]

    #  3- Find the splits that could contain unspliced segments
    #  list of sets. One set per position
    idxs_list = [set([]) for x in range(len(list(records.values())[0].seq) - kmer)]
    aln_length = len(idxs_list)
    if aln_length < min_aln_length:
        exception = "WARNING!: alignment {} is shorter than the allowed \n\
                     minimum STEP3 did not complete".format(
            fasta
        )
        return (exception, "None")
    try:
        singletons_per_index = find_singletons_per_position(
            taxa, aln_length, kmer, sensitivity, alleles
        )
    except:
        exception = "[Error] There is a problem with fasta {}\n\
                     STEP3 did not complete".format(
            fasta
        )
        return (exception, "None")

    split_dict = make_split_dict(singletons_per_index, inverted_taxa_dict)
    problematic_sequences = split_dict.keys()
    corrected_sequences = {}
    l = min_length_fragment

    #  4- Some seqs contain splits that are not shared with other seqs
    #      Investigate whether they are group into stretches longer than
    #      the allowed maximum
    for taxon in problematic_sequences:
        seq = records[taxon].seq
        (start, end) = find_start_end(seq, kmer)
        idxs_to_correct_potential = [
            x for x in sorted(split_dict[taxon]) if x >= start and x <= end
        ]
        idxs_to_correct_intermediate = []
        if idxs_to_correct_potential:
            for (idx, position) in enumerate(idxs_to_correct_potential):
                p = position
                if p <= len(seq) - l:
                    if len(
                        set(idxs_to_correct_potential[idx : idx + l]).intersection(
                            set(range(p, p + l))
                        )
                    ) > differences and "-" * ambiguities not in str(seq[p : p + l]):
                        idxs_to_correct_intermediate.extend(range(p, p + l - kmer))

            if idxs_to_correct_intermediate:
                final_sequence = []
                final_index_to_correct = list(set(idxs_to_correct_intermediate))
                if final_index_to_correct:
                    for idx in range(len(seq)):
                        if idx in final_index_to_correct:
                            final_sequence.append("N")
                        else:
                            final_sequence.append(seq[idx])
                    corrected_sequences[taxon] = "".join(final_sequence)
    not_corrected_sequences = [
        taxon
        for taxon in records.keys()
        if taxon not in list(corrected_sequences.keys())
    ]
    for taxon in not_corrected_sequences:
        corrected_sequences[taxon] = str(records[taxon].seq)

    #  5- Remove sequences containing more ambiguious positions
    #     than the allowed maximum
    remove_taxon = []
    for taxon in corrected_sequences.keys():
        c = Counter(corrected_sequences[taxon])
        ambiguity = c["N"]
        non_ambiguous = sum([c[x] for x in list("aAgGtTcC")])
        if non_ambiguous == 0:
            remove_taxon.append(taxon)
        elif ambiguity / non_ambiguous > max_ambiguity:
            remove_taxon.append(taxon)
    for taxon in remove_taxon:
        del corrected_sequences[taxon]
    if (
        len(list(set([x.split("|")[0] for x in list(corrected_sequences.keys())])))
        >= min_taxa
    ):
        write_fasta(basename + out_suffix + "_complete.fasta", corrected_sequences)
        try:
            shutil.move(basename + out_suffix + "_complete.fasta", out_folder)
        except Exception as ex:
            s = "move fasta: {}".format(fasta)
            message = error_template.format(type(ex).__name__, ex.args, s)
            print(message)
            sys.exit(0)
    else:
        write_fasta(basename + out_suffix + "_incomplete.fasta", corrected_sequences)
        try:
            shutil.move(basename + out_suffix + "_incomplete.fasta", out_folder)
        except Exception as ex:
            s = "move fasta: {}".format(fasta)
            message = error_template.format(type(ex).__name__, ex.args, s)
            print(message)
            sys.exit(0)
    s_i = "STEP3: Done with alignment {}".format(fasta)
    return (s_i, "None")

