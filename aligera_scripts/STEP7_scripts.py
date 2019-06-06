#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Tools that assembled sequences fragments into full size alleles.
    1. Distances (corrected or uncorrected according to user selection) are computed 
        between each sequence pair using the EMBOSS distmat software.
    2. Find all putative alleles. Two sequences are recognized as alleles if they have
        large overlaps, if they are closer together than they are to any other 
        sequence and if their distance is smaller than the user set maximum for allele
        separation.
    3. For taxa that possess several non (low) overlapping sequences that are putative
        alleles, find beginning and end for each sequence, defined as a kmer length 
        region without ambiguity.
    4. Build a directed weighted graph for each taxon with putative alleles, with 
        sequences as edges. An left to right edge is inserted if two sequences have 
        little or no overlap and if they share the closest overlapping sequence from a 
        different taxon. 
    5. Brute force algorithm that generate all possible path between nodes that have 
        outdegrees of one and node that have indegree of one. Path are scored
        according to their length in term of number of nucleotids minus the number of 
        weighted mismaches, minus the number of weighted overlap.
    6. The path with the highest score is is selected and nodes involved in the path
        are removed. Path are generated from the remaining nodes and the best path is
        selected. These operations are repeated  until all nodes have been used at 
        least once.
"""
__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import os
import shutil
import networkx as nx
from itertools import chain, combinations, product
import collections
from math import ceil

from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

try:
    from aligera_scripts.utilities import (
        sequence_profiler,
        get_overlap_score,
        calculate_distance_matrix,
        distance_matrix_parser,
    )
except:
    from utilities import (
        sequence_profiler,
        get_overlap_score,
        calculate_distance_matrix,
        distance_matrix_parser,
    )


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def find_start_end(seq, kmer):
    """
    Finds the start and end of a sequence in an alignment,
    identified by kmer consecutive non "-" bases.
    """
    start, end = None, None
    for idx in range(len(seq)):
        if seq[idx] != "-" and "-" not in seq[idx : idx + kmer]:
            start = idx
            break
    for idx in range(len(seq)):
        if seq[-idx] != "-" and "-" not in reversed(seq[-idx : -idx - kmer]):
            end = len(seq) - idx
            break
    if start < end:
        return start, end


def get_distance(item, dist_dict):
    """
    Return an infiny value if the key is not in the dictionary.
    Function for closest_sequence_match() and assertain_allele().
    """
    try:
        d = dist_dict[item]
    except:
        d = float("inf")
    return d


def assertain_allele(
    taxa_contigs,
    dist_dict,
    taxa_profiles,
    dict_overlap_seq,
    max_alleles_dist,
    overlapping_seq,
):
    """
    Function that evaluate if two sequence can be alleles
    It check that two sequences are closer together than they are to any other sequence.
    """
    pairs = [tuple(sorted(x)) for x in combinations(taxa_contigs, 2)]
    overlap_pairs = [p for p in pairs if p[0] in dict_overlap_seq[p[1]]]
    association_dict = {}

    for p in overlap_pairs:
        most_similar_0 = sorted(
            [
                (item, dist_dict[tuple(sorted([p[0], item]))])
                for item in dict_overlap_seq[p[0]]
                if (
                    get_distance(tuple(sorted([p[0], item])), dist_dict)
                    <= max_alleles_dist
                )
            ],
            key=lambda x: x[1],
        )

        for match in most_similar_0[:overlapping_seq]:
            if match[0] == p[1]:
                d = match[1]
                if p[0] not in list(association_dict.keys()) and p[1] not in list(
                    association_dict.keys()
                ):
                    association_dict[p[0]] = [p[1], d]
                    association_dict[p[1]] = [p[0], d]

                elif p[0] not in list(association_dict.keys()) and p[1] in list(
                    association_dict.keys()
                ):
                    if d < association_dict[p[1]][1]:
                        previous_association_1 = association_dict[p[1]][0]
                        association_dict[p[0]] = [p[1], d]
                        association_dict[p[1]] = [p[0], d]
                        del association_dict[previous_association_1]

                elif p[1] not in list(association_dict.keys()) and p[0] in list(
                    association_dict.keys()
                ):
                    if d < association_dict[p[0]][1]:
                        previous_association_0 = association_dict[p[0]][0]
                        association_dict[p[0]] = [p[1], d]
                        association_dict[p[1]] = [p[0], d]
                        del association_dict[previous_association_0]

                else:
                    if d < association_dict[p[1]][1] and d < association_dict[p[0]][1]:
                        previous_association_0 = association_dict[p[0]][0]
                        previous_association_1 = association_dict[p[1]][0]
                        association_dict[p[0]] = [p[1], d]
                        association_dict[p[1]] = [p[0], d]
                        del association_dict[previous_association_0]
                        del association_dict[previous_association_1]
                break

    #  Keep a the longest alleles in the allele pair
    final_association = []
    for p in list(set([tuple(sorted([k, v[0]])) for k, v in association_dict.items()])):
        if taxa_profiles[p[0]] > taxa_profiles[p[1]]:
            final_association.append((p[0], p))
        else:
            final_association.append((p[1], p))

    final_association.extend(
        [(x, None) for x in taxa_contigs if x not in list(association_dict.keys())]
    )

    #  The pairs that are not recognized as alleles are deemed to be paralogs
    para_dict = {}
    associated_alleles = {x: sorted(y) for (x, y) in final_association if y is not None}

    for p in pairs:
        s_p = sorted(p)
        para_keys = list(para_dict.keys())
        assoc_alleles = list(associated_alleles.keys())
        if s_p not in list(associated_alleles.values()):
            if p[0] not in para_keys and p[0] in assoc_alleles:
                para_dict[p[0]] = [p[1]]
            elif p[0] in para_keys and p[0] in assoc_alleles:
                para_dict[p[0]].extend([p[1]])

            if p[1] not in para_keys and p[1] in assoc_alleles:
                para_dict[p[1]] = [p[0]]
            elif p[1] in para_keys and p[1] in assoc_alleles:
                para_dict[p[1]].extend([p[0]])
    return final_association, para_dict


def best_matches(distances_to_pair, overlap_ratio):
    """
    Find the closest seq from  a target sequence (in a different taxon),
    allele pairs are counted as a single sequence.
    overlap_ratio corresponds to the nbr of accepted_seq divided
    by the total number of sequences that overlap the candidate pair.
    Function for closest_sequence_match().
    """
    #  number of seq from other taxa that overlap the target seq
    total_overlap = len(set([x[0].split("allele")[0] for x in distances_to_pair]))
    #  index in list needed to reach the overlap_ratio, float rounded to next int
    list_idx = ceil(total_overlap * overlap_ratio)
    #  return the list of closest seq (without the distance) of length list_idx
    return [x[0] for x in distances_to_pair[:list_idx]]


def closest_match(
    p,
    taxon,
    align_dict,
    taxa_profiles,
    dist_dict,
    overlap_threshold,
    overlap_ratio,
    alleles_dict,
):
    """
    Finds overlapping sequences of both members in the pair and evaluate
    whether these sequences are in the common set of closest sequences.
    Function for alleles_assignment().
    """
    #  Alleles pair should only be represented by the longest allele.
    #  alleles_dict_rev reverses the k, v relationship of alleles_dict
    #  so each allele is a k and the longest allele of the pair is the v.
    #  1- add the seqs that are not allele
    alleles_dict_rev = {k: k for k, v in alleles_dict.items() if v is not None}
    #  2- add the seqs that are alleles
    for k, v in alleles_dict.items():
        if v is not None:
            for item in v:
                alleles_dict_rev[item] = k

    #  The quandidate seq take into account only the longest of the two alleles
    quandidate_seq = set(
        [
            alleles_dict_rev[x]
            for x in align_dict.keys()
            if x not in p
            and x.split("|")[0] != taxon
            and x in list(alleles_dict_rev.keys())
        ]
    )
    #  Seqs from other taxa that overlap both seqs in pair
    overlap_seqs = [
        x
        for x in quandidate_seq
        if (get_overlap_score(sorted([x, p[0]]), taxa_profiles) > overlap_threshold / 2)
        and (
            get_overlap_score(sorted([x, p[1]]), taxa_profiles) > overlap_threshold / 2
        )
    ]
    if overlap_seqs:
        dist_p_0 = sorted(
            [
                (x, get_distance(tuple(sorted([p[0], x])), dist_dict))
                for x in overlap_seqs
            ],
            key=lambda x: x[1],
        )

        best_matches_0 = best_matches(dist_p_0, overlap_ratio)

        dist_p_1 = sorted(
            [
                (x, get_distance(tuple(sorted([p[1], x])), dist_dict))
                for x in overlap_seqs
            ],
            key=lambda x: x[1],
        )

        best_matches_1 = best_matches(dist_p_1, overlap_ratio)
        if set(best_matches_0).intersection(set(best_matches_1)):
            return True
    return False


def make_consensus(path, align):
    """
    Concatenate sequences in path.
    Differences are resolved by taking the nucleotide in majority within the column.
    """
    #  Alignment of sequences in the path
    path_aln = MultipleSeqAlignment([rec for rec in align if rec.name in path])

    # Alignment of sequences NOT in the path
    no_path_aln = MultipleSeqAlignment([rec for rec in align if rec.name not in path])

    consensus_sequence = ""
    for idx in range(path_aln.get_alignment_length()):
        path_col = path_aln[:, idx]
        ambiguity = set([x for x in path_col if x != "-"])
        if len(ambiguity) > 1:
            no_path_col = [x for x in no_path_aln[:, idx] if x != "-"]
            if no_path_col:
                c = collections.Counter(sorted(no_path_col))
                fq = sorted(
                    [(x, c[x] / len(no_path_col)) for x in c.keys()],
                    key=lambda x: x[1],
                    reverse=True,
                )

                shared_major_base = set([x[0] for x in fq if x[1] > 0.25]).intersection(
                    ambiguity
                )
                ambiguity_up = set([x.upper() for x in ambiguity])

                #  Case 1 one of the two variants is among the dominant bases,
                #    ambiguity is resolved by using the dominant base
                if shared_major_base and len(shared_major_base) == 1:
                    consensus_sequence += list(shared_major_base)[0]

                #  Case 2 none of the variants are among the dominant bases or
                #   both of them are among the dominant bases,
                #   then return the consensus
                elif ambiguity_up in [
                    set(["A", "G"]),
                    set(["C", "T"]),
                    set(["G", "C"]),
                    set(["A", "T"]),
                    set(["G", "T"]),
                    set(["A", "C"]),
                    set(["C", "G", "T"]),
                    set(["A", "G", "T"]),
                    set(["A", "C", "T"]),
                    set(["A", "C", "G"]),
                    set(["A", "T", "C", "G"]),
                ]:
                    iupac_code = {
                        ("A", "G"): "R",
                        ("C", "T"): "Y",
                        ("G", "C"): "S",
                        ("A", "T"): "W",
                        ("G", "T"): "K",
                        ("A", "C"): "M",
                        ("C", "G", "T"): "B",
                        ("A", "G", "T"): "D",
                        ("A", "C", "T"): "H",
                        ("A", "C", "G"): "V",
                        ("A", "T", "C", "G"): "N",
                    }
                    iupac_ambiguity = [
                        iupac_code[x]
                        for x in iupac_code.keys()
                        if set(x).intersection(ambiguity_up) == set(x)
                    ][0]
                    consensus_sequence += iupac_ambiguity

                # Case 3 the alleles already contain ambiguities
                else:
                    consensus_sequence += "N"
        else:
            if [x for x in path_col if x != "-"]:
                consensus_sequence += [x for x in path_col if x != "-"][0].replace(
                    "N", "-"
                )
            else:
                consensus_sequence += "-"
    new_name = path[0].split("|")[0] + "|" + "|".join([x.split("|")[1] for x in path])

    return (new_name, consensus_sequence)


def get_path_length(path, start_end_dict):
    """
    Calculate the length of alleles (nbr of bases) along a path
    """
    seq_limits = [start_end_dict[x[0]] for x in path]
    S = set([])
    for pair in seq_limits:
        S = S.union(set(range(pair[0], pair[1])))
    return len(S)


def count_mismatches_overlap_old(path, records_all):
    """
    Calculate the amout of mismatches along a path
    """
    records_with_alleles = records_all.takeSeqs([x[0] for x in path])
    i = 0
    j = 0
    for column in records_with_alleles.iterPositions():
        M = [x for x in column if x != "-"]
        N = set(M)
        if len(M) > 1:
            i += 1
        if len(N) > 1:
            j += 1
    return i, j


def count_mismatches(path, align):
    """
    Calculate the amout of mismatches along a path
    """
    names = [x[0] for x in path]
    sliced_align = MultipleSeqAlignment([rec for rec in align if rec.name in names])
    i = 0
    j = 0
    for idx in range(sliced_align.get_alignment_length()):
        column = sliced_align[:, idx]
        M = [x for x in column if x != "-"]
        N = set(M)
        if len(M) > 1:
            i += 1
        if len(N) > 1:
            j += 1
    return i, j


def rename(seq_name):
    """
    Rename sequences. 'contig' is removed and 'allele' is replace by 'a' 
    """
    new_name = seq_name.replace("allele", "a").replace("contig_", "")
    return new_name


def assemble_alleles(fasta, cfg, **kargs):
    """
    Main function for STEP7:
    Allele fragments are assembled in order to recover the full length transcript.
    interp. 
        overlap_threshold (int)         Parameter that plays two roles in the analysis:
                                        1. in order to investigate two fragment from 
                                        the same origin for they cannot overlap by more 
                                        than this value.
                                        2. in order to assemble two fragment f0, f1, 
                                        there must be a fragment t from a different 
                                        taxon that overlap each of the two sequences by 
                                        overlapping_threshold / 2. f0 and f1 are 
                                        assembled if they are both in the overlap region
                                        among the closest sequences from t.
        kmer (int)                      The number of non ambiguous consecutive bases to 
                                        mark the beginning and end of a sequence. 
        overlap_proportion (float)      When testing for homology of two non (low) 
                                        overlapping sequences f0, f1 from taxon t0; 
                                        we use a search the set of sequences from the
                                        other taxa that largely overlap both f0 and f1.
                                        The number of taxa to be search is 
                                        (number seq that overlap) * overlap_proportion
        max_alleles_distance (float)    indicates the max distance allowed between 2 
                                        alleles. 
        mistmatch_cost (int)            Penalty cost for mismatch in the overlapping 
                                        flanking regions of putative alleles.
        overlap_cost (int)              Penalty cost for overlap in the flanking regions 
                                        of putative alleles.                                    
    
    """
    in_format = cfg["input_format"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    taxa_without_alleles = cfg["taxa_without_alleles"]
    if taxa_without_alleles:
        taxa_no_allele = taxa_without_alleles
    else:
        taxa_no_allele = []
    kmer = cfg["kmer"]
    overlap_ratio = cfg["overlap_proportion"]
    overlap_threshold = cfg["overlap_threshold"]
    max_allele_dist = cfg["max_alleles_distance"]
    overlap_cost = cfg["overlap_cost"]
    mistmatch_cost = cfg["mistmatch_cost"]
    basename = fasta.split(in_format)[0]

    try:
        assert (kmer) < 50
    except:
        s = "'kmer' > 50, unrealistic parameter value"
        raise ValueError(s)
        sys.exit(s)

    try:
        cfg["dna_model"] in [0, 1, 2, 3, 4, 5]
    except:
        s = "Wrong parameters: 'dna_model' must be in [0, 1, 2, 3, 4, 5]"
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

    align_length = len(align)
    if align_length < 3:
        s_i = "STEP7 Done with fasta: {}".format(fasta)
        s_d = "No transformation performed on fasta {0} \n\
which only contains {1} sequences".format(
            fasta, align_length
        )
        return (s_i, s_d)

    length = align.get_alignment_length()
    try:
        assert overlap_threshold < length
    except:
        exception = "[Error]: for fasta {} 'overlap_threshold'\
         exceeds the length of the alignment".format(
            fasta
        )
        return (exception, None)

    #  Calculate distance matrix
    calculate_distance_matrix(
        fasta, {"input_suffix": in_format, "dna_model": cfg["dna_model"]}
    )
    distance_matrix = basename + "_for_dnadist.fasta.distmat"
    dist = distance_matrix_parser(distance_matrix)
    dist_dict = {tuple(sorted([x[0], x[1]])): x[2]["distance"] for x in dist}
    non_empty_seqs = list(set(list(chain(*list(dist_dict.keys())))))

    #  Find records with alleles that have not been excluded in the settings
    rec_alleles = [
        rec
        for rec in align
        if rec.name.split("|")[0] not in taxa_no_allele and rec.name in non_empty_seqs
    ]
    align_dict = {rec.name: str(rec.seq) for rec in align if rec.name in non_empty_seqs}
    taxa = {rec.name: str(rec.seq) for rec in rec_alleles if rec.name in non_empty_seqs}

    #  For each sequence we recode 0 for '-', 1 otherwise,
    #   into a dictionary {name_number: profile}"
    taxa_profiles = sequence_profiler(taxa)

    #  Find beginning and end for each sequence.
    #  If no/beginning end is found do not include the sequence
    start_end_dict = {
        rec: (find_start_end(taxa[rec], kmer))
        for rec in taxa
        if None not in [find_start_end(taxa[rec], kmer)]
    }

    #  Nbr ovelapping sequences used for calculation
    overlapping_seq = int(overlap_ratio * len(list(taxa.keys())))

    #  Dictionary that contains the taxa (i.e. the transcriptomes) as keys
    #   and the names of all sequences corresponding to the key as value.
    taxa_list = list(set([name.split("|")[0] for name in taxa.keys()]))
    taxa_contigs = dict(zip(taxa_list, [[] for x in taxa_list]))
    for rec in taxa:
        if rec in list(start_end_dict.keys()):
            taxon = rec.split("|")[0]
            taxa_contigs[taxon].append(rec)

    #  Dict with items of overlapping_seq as keys and as values the list
    #   of all overlapping seqs. Dict used for finding the seqs from
    #   a different taxon  that overlap the sequences from a given taxon.
    dict_overlap_seq = {}
    for seq in taxa:
        dict_overlap_seq[seq] = [
            x
            for x in taxa.keys()
            if x != seq
            and (
                get_overlap_score(tuple(sorted([x, seq])), taxa_profiles)
                > overlap_threshold
            )
        ]
    #  Obtain a dictionary that keep for each allele pair the longest and for
    #   the sequences that are not paired assigns the key 'None'.
    alleles_dict = {}
    para_dict = {}
    for taxon in [x.split("|")[0] for x in taxa.keys()]:
        if len(taxa_contigs[taxon]) > 1:
            L, P = assertain_allele(
                taxa_contigs[taxon],
                dist_dict,
                taxa_profiles,
                dict_overlap_seq,
                max_allele_dist,
                overlapping_seq,
            )
            para_dict.update(P)
            for (k, v) in L:
                alleles_dict[k] = v
        else:
            alleles_dict[taxa_contigs[taxon][0]] = None
    #  Obtain a new dictionary of taxon as key and sequences without alleles
    #   as values. These are the sequences that are not alleles so they could
    #   be assembled into full length transcripts.
    taxa_contigs_no_alleles = {}
    for k, v in list(taxa_contigs.items()):
        taxa_contigs_no_alleles[k] = [
            seq for seq in v if seq in list(alleles_dict.keys())
        ]

    #  Dict that for each taxon find the seq pairs with overlap below the
    #   threshold and that have a common sequences that overlap both sequences
    #   as the best match as obtained from the distance analysis.
    taxa_with_pair = {}
    out_seqs = []
    used_seq_in_out = []

    for taxon, seqs in list(taxa_contigs_no_alleles.items()):

        #  Case 1. there is only one sequence, so no pair can be infered.
        if len(seqs) == 1:
            #  Case the unique sequence has an allele pair
            if alleles_dict[seqs[0]] and len(alleles_dict[seqs[0]]) > 1:
                out_seqs.extend(
                    [
                        (x, align_dict[x].upper().replace("N", "-"))
                        for x in alleles_dict[seqs[0]]
                    ]
                )
                used_seq_in_out.extend([x for x in alleles_dict[seqs[0]]])
            #  Case the unique sequence has no allele pair
            else:
                out_seqs.append(
                    (seqs[0], align_dict[seqs[0]].upper().replace("N", "-"))
                )
                used_seq_in_out.append(seqs[0])
            continue

        #  Case 2. there are several sequences, try to find ways to concatenate them.
        seqs_pairs = list(combinations(seqs, 2))
        L = []
        for p in seqs_pairs:
            overlap = get_overlap_score(p, taxa_profiles)
            start_0, end_0 = start_end_dict[p[0]]
            start_1, end_1 = start_end_dict[p[1]]
            s_p = sorted([(p[0], start_0), (p[1], start_1)], key=lambda x: x[1])
            sorted_pair = (s_p[0], s_p[1])
            if not overlap:
                seq_coverage = (end_0 - start_0) + (end_1 - start_1)
            else:
                seq_coverage = (
                    max(end_0, end_1)
                    - min(end_0, end_1)
                    + min(end_0, end_1)
                    - max(start_0, start_1)
                    + max(start_0, start_1)
                    - min(start_0, start_1)
                )
            #  If they are little overlap between two seqs of the same taxon
            #   and they share the closest overlapping sequence, create
            #   an edge between the two sequences.
            if overlap <= overlap_threshold and closest_match(
                p,
                taxon,
                align_dict,
                taxa_profiles,
                dist_dict,
                overlap_threshold,
                overlap_ratio,
                alleles_dict,
            ):
                L.append(
                    (
                        sorted_pair,
                        {
                            "overlap": overlap,
                            "sequence_coverage": seq_coverage,
                            "start": min(start_0, start_1),
                        },
                    )
                )
        if L:
            taxa_with_pair[taxon] = sorted(L, key=lambda x: x[1]["start"])
        else:
            #  No pair could be calculated
            for p in seqs_pairs:
                for item in p:
                    if alleles_dict[item] and len(alleles_dict[item]) > 1:
                        out_seqs.extend(
                            [
                                (x, align_dict[x].upper().replace("N", "-"))
                                for x in alleles_dict[item]
                            ]
                        )
                        used_seq_in_out.extend([x for x in alleles_dict[item]])
                    else:
                        out_seqs.append(
                            (item, align_dict[item].upper().replace("N", "-"))
                        )
                        used_seq_in_out.append(alleles_dict[item])
                continue

    #  Case there are no sequences that can be paired:
    #  Return the original alignment.
    if not taxa_with_pair:
        clean_align = [rec for rec in align if rec.name in non_empty_seqs]
        out_name = basename + out_suffix + ".fasta"
        SeqIO.write(clean_align, out_name, "fasta")
        if not os.path.isfile(os.path.join(out_folder, out_name)):
            shutil.move(out_name, out_folder)
        s_i = "STEP7 Done with fasta: {}".format(fasta)
        s_d = "No full size allele could be reconstructed for fasta: {}".format(fasta)

        return (s_i, s_d)

    #  Case some sequences are paired:
    #  1. Generate all possible path between nodes that have indegrees of one
    #     and nodes with with outdegree of one.
    #  2. The path are scored according to their length in term of number
    #     of nucleotids minus the number of weighted mismaches, minus
    #     the number of weighted overlap.
    #  3. The best path is selected.
    #  4. Then new path are recalculated by removing the selected terminal
    #     nodes and recalculating the path. Adding then to the path pool
    #     for selection.
    #  5. Continue until all nodes have been used at least once.
    for taxon in list(taxa_with_pair.keys())[:]:
        G = nx.DiGraph()
        G.add_edges_from([(x[0][0], x[0][1], x[1]) for x in taxa_with_pair[taxon]])
        in_nodes = [k for k, v in list(G.in_degree(G.nodes())) if v == 0]
        out_nodes = [k for k, v in list(G.out_degree(G.nodes())) if v == 0]

        #   Initialize the search procedure with the first path
        possible_paths = []
        for item in product(in_nodes, out_nodes):
            possible_paths.extend(
                nx.all_simple_paths(G, source=item[0], target=item[1])
            )
        used_seqs = []
        possible_path_dict = {}
        for path in possible_paths:
            score = (
                get_path_length(path, start_end_dict)
                - overlap_cost * count_mismatches(path, align)[0]
                - mistmatch_cost * count_mismatches(path, align)[1]
            )
            possible_path_dict[tuple(path)] = score

        best_path = sorted(list(possible_path_dict.items()), key=lambda x: x[1])[-1]
        sequences_best_path = [y for y in best_path[0]]

        all_best_paths = [best_path]
        remaining_seqs = [x for x in G.nodes() if x not in [y[0] for y in best_path[0]]]

        #  Remove the sequences of the best path from the graph only
        #   if there is a paralog sequence that is overlapping them.
        for seq in sequences_best_path:
            used_seqs.append(seq[0])
            if seq[0] in list(para_dict.keys()):
                G.remove_node(seq)
                for path in list(possible_path_dict.keys()):
                    if seq in path:
                        del possible_path_dict[path]

        remaining_seqs = list(G.nodes())

        while len(remaining_seqs) > 1:
            #  Retrieve 5' and 3' sequences.
            in_nodes = [k for k, v in list(G.in_degree(G.nodes())) if v == 0]
            out_nodes = [k for k, v in list(G.out_degree(G.nodes())) if v == 0]
            if in_nodes == out_nodes:
                break

            #  Find all path that originate in one in_node and end in one out_node
            possible_paths = []
            for item in product(in_nodes, out_nodes):
                possible_paths.extend(
                    nx.all_simple_paths(G, source=item[0], target=item[1])
                )
            for path in possible_paths:
                score = (
                    get_path_length(path, start_end_dict)
                    - overlap_cost * count_mismatches(path, align)[0]
                    - mistmatch_cost * count_mismatches(path, align)[1]
                )
                possible_path_dict[tuple(path)] = score

            #  Select the path with the best score:
            #   if the putative path is entirely included in a previously
            #   selected path, it is rejected and the loop ends
            for path in sorted(
                list(possible_path_dict.items()), key=lambda x: x[1], reverse=True
            ):
                for b_path in [set(x[0]) for x in all_best_paths]:
                    if set(path[0]).intersection(b_path) == set(path[0]):
                        best_path = None
                        break
                if path not in all_best_paths and best_path is not None:
                    best_path = path
                    break
                else:
                    best_path = None
            if best_path is None:
                break

            all_best_paths.append(best_path)
            sequences_best_path = [y for y in best_path[0]]
            #  Switch that records whether the path is entirely made of
            #   sequences lacking paralogs
            no_paralogue = True

            #  1- Remove all sequences in the best path that have a paralog
            for seq in sequences_best_path:
                used_seqs.append(seq[0])
                if seq[0] in list(para_dict.keys()):
                    G.remove_node(seq)
                    used_seqs.append(seq)
                    no_paralogue = False
                    for path in list(possible_path_dict.keys()):
                        if seq in path:
                            del possible_path_dict[path]

            #  2- Remove all sequences in the path if only made out of
            #     sequences lacking a paralog
            if no_paralogue:
                for seq in sequences_best_path:
                    G.remove_node(seq)
                    used_seqs.append(seq)
                    for path in list(possible_path_dict.keys()):
                        if seq in path:
                            del possible_path_dict[path]
            remaining_seqs = G.nodes()

        if remaining_seqs and remaining_seqs[0][0] not in used_seqs:
            all_best_paths.extend([((remaining_seqs[0],), None)])

        all_best_paths_seq = []
        for x in all_best_paths:
            all_best_paths_seq.append([y[0] for y in x[0]])

        final_paths = []
        for path in all_best_paths_seq:
            final_paths.append(path)
            alleles_path = []
            #  Switch that ensures that there is at least one seq with alleles
            switch = False
            for p in path:
                if alleles_dict[p] is None and len(path) > 1:
                    alleles_path.append(p)
                elif alleles_dict[p] is not None:
                    alternative_allele = [x for x in alleles_dict[p] if x != p]
                    alleles_path.extend(alternative_allele)
                    switch = True
            if switch:
                final_paths.append(alleles_path)

        for path in final_paths:
            used_seq_in_out.extend(path)
            out_seqs.append(make_consensus(path, align))

    #  Add sequences that have not been included in any path

    for item in list(
        set([x for x in align_dict.keys()]).difference(set(used_seq_in_out))
    ):
        if item not in [x[0] for x in out_seqs]:
            out_seqs.append((item, align_dict[item].replace("N", "-")))
    out_seq_newname = zip(
        [rename(x[0]) for x in sorted(list(set(out_seqs)), key=lambda x: x[0])],
        [x[1] for x in sorted(list(set(out_seqs)), key=lambda x: x[0])],
    )
    final_records = []
    for name, seq in out_seq_newname:
        final_records.append(
            SeqRecord(
                Seq(seq, IUPAC.IUPACAmbiguousDNA()), name=name, id=name, description=""
            )
        )

    out_name = basename + out_suffix + ".fasta"

    SeqIO.write(final_records, out_name, "fasta")
    if not os.path.isfile(os.path.join(out_folder, out_name)):
        shutil.move(out_name, out_folder)

    s_i = "STEP7 Done with fasta: {}".format(fasta)
    s_d = "Alignment size was reduced from {0} to {1} sequences".format(
        len(align), len(final_records)
    )
    return (s_i, s_d)


