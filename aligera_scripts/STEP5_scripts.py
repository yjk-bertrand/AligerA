#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Tools that separates outparalogs.
    1. Remove empty sequences.
    2. Alleles that have been phased are identified by their 'allele' tag name. 
        When several sequences from the same taxon overlap and when 'alleles' have been
        phased, we identify "alleles" and putative paralogs.
    3. Distances (corrected or uncorrected according to user selection) are computed 
        between each sequence pair using the EMBOSS distmat software.
    4. For each taxon, distances between putative paralogs are computed. The minimum
        of these distances (d) is retained.
    5. A weighted undirected graph is built with sequences as vertices and distance as 
        edges.
        Edges longer than than d are removed yielding connected commponents. 
        If a single component is obtained, return the original alignment.
    6. When several components have been obtained, we test for overseparation by 
        computing the the Jaccard index for the number of taxa in the components:
        J(C1,C2) = (taxa in C1 and C2)/(taxa in C1 or C2). If the index's value
        is below a present threshold the sequence data is check to see if fusion of 
        components does not introduce paralogs.
    7. If there is more than two components or two components and ungroup sequences,
        the smallest components (less than half the size of the largest) are 
        desassembled. Their sequences and the ungroup sequences are assigned to the 
        larger components based on smaller distance to one of the sequences in the 
        components distance as long as the distance are below the maximum distance
        threshold set by the user.
    8. Files that harbour less taxa than the set threshold are tagged as '_incomplete',
        otherwise they are tagged as '_complete'. Sequences that remained unassigned
        after the distance clustering are labelled as 'unassigned_seq'..
    
"""
__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import os
import shutil
import operator
import itertools

from Bio import SeqIO

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

try:
    from aligera_scripts.STEP4_scripts import (
        does_overlap,
        calculate_Jaccard_index,
        write_fasta,
    )
except:
    from STEP4_scripts import does_overlap, calculate_Jaccard_index, write_fasta

try:
    from aligera_scripts.STEP7_scripts import get_distance
except:
    from STEP7_scripts import get_distance


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def dist_assignment(
    components,
    records,
    taxa_profiles,
    distances,
    basename,
    out_folder,
    max_dist,
    min_taxa,
    outgroups,
    min_overlap,
    out_suffix,
):
    """
    Sequences that do not cluster into components are assigned based on distance method
    This version contains small differences with the STEP4 function:
        it does not separate small and large components
    """
    comps = components
    final_comps = {tuple(comp): comp for comp in comps}
    unassigned_seqs = list(set(records.keys()).difference(set(itertools.chain(*comps))))
    #
    #  If there is at least one unassigned sequence
    #
    if len(unassigned_seqs) > 0:
        #  Use distances to assign remaining taxa'
        if not distances:
            exception = "[Error] problem parsing distances\nExiting"
            raise Exception(exception)
            sys.exit()
        unassignable_taxa = []
        for unassigned_seq in unassigned_seqs:
            min_dist = {}
            for comp in comps:
                D = sorted(
                    [
                        x
                        for x in distances
                        if (unassigned_seq in x[0] and x[1] in comp)
                        or (unassigned_seq in x[1] and x[0] in comp)
                        and x[2]["distance"] != "nan"
                    ],
                    key=lambda x: x[2]["distance"],
                )

                if D:
                    min_dist[tuple(comp)] = D[0][2]["distance"]
            sorted_min_dist = sorted(min_dist.items(), key=operator.itemgetter(1))

            #  If one of the components already contains the unassigned taxon,
            #  only a limited overlap between the two sequences is allowed
            if sorted_min_dist:

                #  List of pairs containing one unassigned seq and
                #  one seq from the selected component
                overlaps = [[unassigned_seq, x] for x in sorted_min_dist[0][0]]

                #   Calculate overlap between unassigned seq and each member
                #   of component. Return the max.
                max_overlap = max(
                    [
                        sum(
                            a == b
                            for a, b in zip(
                                taxa_profiles[item[0]], taxa_profiles[item[1]]
                            )
                            if a == 1
                        )
                        for item in overlaps
                    ]
                )

                #  Case 1. There is overlap with seq from the same taxon and
                #  the dist is smaller than the threshold
                #  or the overlap is small but the distance is minute
                if unassigned_seq.split("|")[0] not in outgroups:
                    if (
                        sorted_min_dist[0][1] <= max_dist and max_overlap > min_overlap
                    ) or sorted_min_dist[0][1] <= max_dist / 3:
                        #  Unless there is already a sequence that overlapp
                        new_final_comps = final_comps[sorted_min_dist[0][0]][:]
                        new_final_comps.append(unassigned_seq)
                        final_comps[sorted_min_dist[0][0]] = new_final_comps
                    else:
                        unassignable_taxa.append(unassigned_seq)

                #  Case 2. We are dealing with the outgroup and
                #  we double the max dist for assigning seq to component
                else:
                    if (
                        sorted_min_dist[0][1] <= max_dist * 2
                        and max_overlap > min_overlap
                    ) or sorted_min_dist[0][1] <= max_dist / 3:
                        new_final_comps = final_comps[sorted_min_dist[0][0]][:]
                        new_final_comps.append(unassigned_seq)
                        final_comps[sorted_min_dist[0][0]] = new_final_comps
                    else:
                        unassignable_taxa.append(unassigned_seq)

        if unassignable_taxa:
            unassignable_comps = {tuple(unassignable_taxa): unassignable_taxa}
            write_fasta(
                unassignable_comps,
                records,
                basename,
                out_folder,
                min_taxa,
                out_suffix + "_unassigned_seq",
            )
    return final_comps


def get_dist_other_taxa(seq, dist_dict, max_dist, dict_overlap_seq):
    """
    for a seq produces a list of dists. between the seq. and 
    overlapping seqs. from other taxa if the distance (dist) is below the max_dist       
    interp.: 
            seq              (string)    the name of the target seq. 
            dist_dict        (dict)      holds the dists. Pair of seqs. as keys, 
                                         dist. as values.
            max_dist         (int)       maximum allowed dist between 2 seqs. 
            dict_overlap_seq (dict)      holds seqs. as keys and  as values list of 
                                         overlapping seqs from a different taxa.                                    
    """

    return [
        dist_dict[tuple(sorted([seq, p]))]
        for p in dict_overlap_seq[seq]
        if (get_distance(tuple(sorted([seq, p])), dist_dict) < max_dist)
    ]


def inpara_validation(
    taxa,
    taxa_profiles,
    min_overlap,
    taxa_no_outgp_para,
    dists,
    dist_ratio,
    dict_overlap_seq,
    max_alleles_dist,
):
    """
    Function that check whether the alignment, genuinely contains
    a sufficient amount of sequences grouped by taxa with large
    overlap and large distances. These sequences are then deemed
    to be genuine paralogous sequences.
    interpr. :
    taxa (dict)                 contains taxa that are not outgroups as key
                                  and their seqs. as values.
    taxa_profiles (dict)        contains for each seq its bit array, 
                                  (1 = nucleotide, 0 = "-" or "n" ).
    min_overlap (int)           length of min overlap between 2 seq. required in order 
                                  to use their dist. in comparisons.
    taxa_no_outgp_para (dict)   contains seq. names that neither outgroups .
                                  not excluded taxa as keys and their seqs. as values                             
    dists (dict)                keeps for all seq. pair its computed distance. 
    dict_overlap_seq (dict)     keeps for each seq the list of seqs. from  other taxa
                                  that overlap with it.
    max_alleles_dist (float)    indicates the max distance allowed between 2 alleles . 
                                                       
    """
    n_taxa_names = len(set([x.split("|")[0] for x in list(taxa_no_outgp_para.keys())]))
    split_taxa = sorted(
        [tuple(x.split("|")) for x in list(taxa_no_outgp_para.keys())],
        key=lambda x: x[0],
    )
    taxa_gps = []
    for key, gp in itertools.groupby(split_taxa, lambda x: x[0]):
        taxa_gps.append(list(gp))

    #  Taxa with more than a single sequence
    gps = [x for x in taxa_gps if len(x) > 1]
    #  Turn the distance list, dists, into a distance dictionary
    dist_dict = {tuple(sorted([x[0], x[1]])): x[2]["distance"] for x in dists}
    #  Distances between paralogs
    dist_para = []
    #  Distances between alleles
    dist_alleles = []

    #  Identify the presence of overlapping sequences from the same taxon
    for taxon in gps:
        taxa_overlap = []
        #  Test all sequences combinations from the same taxon
        for pair in itertools.combinations(taxon, 2):
            item = tuple(
                sorted([pair[0][0] + "|" + pair[0][1], pair[1][0] + "|" + pair[1][1]])
            )
            if (
                item in list(dist_dict.keys())
                and does_overlap(item, min_overlap, taxa_profiles)
                and dist_dict[item] != "nan"
            ):
                dist_pair = dist_dict[item]
                if dist_pair > max_alleles_dist:
                    dist_para.append((item[0], item[1], dist_dict[item]))
                else:
                    a0 = item[0]
                    dist_a0_para = get_dist_other_taxa(
                        a0, dist_dict, max_alleles_dist, dict_overlap_seq
                    )
                    a1 = item[1]
                    dist_a1_para = get_dist_other_taxa(
                        a1, dist_dict, max_alleles_dist, dict_overlap_seq
                    )
                    #  Situation when there are no seq from another taxon
                    #  overlapping with the allele pair
                    if not dist_a0_para + dist_a1_para:
                        continue
                    min_dist_para = min(dist_a0_para + dist_a1_para)
                    if dist_pair > min_dist_para:
                        dist_para.append((a0, a1, dist_dict[item]))
                    else:
                        dist_alleles.append((a0, a1, dist_dict[item]))

    #  Case when no putative paralog has been identified
    if not dist_para:
        return None, None, None, None

    sort_dist_para = [x[2] for x in dist_para]
    #  Max dist between two paralogs
    max_dist = max(sort_dist_para)

    #  Retrieve max and min distances for paralog separation
    #  Find within the allele list a taxon that contains alleles in a pair
    #  in the paralogs list. Such allele pair of alleles is considered
    #  to be a genuine pair of alleles if the distance between the alleles is
    #  at least N (user defined) times lower than the distance between
    #  the the paralogs.
    #  The largest of the allele distances is used as a threshold for
    #  separating paralogs from alleles.

    min_dist = 0
    if dist_alleles:
        for pair in dist_alleles:
            hits = []
            for item in dist_para:
                if set([item[0], item[1]]).intersection([pair[0], pair[1]]):
                    hits.append(item[2])
            if hits:
                paralogs_max_dist = max(hits)
                if paralogs_max_dist > pair[2] * dist_ratio and pair[2] > min_dist:
                    min_dist = pair[2]
    else:
        min_dist = sort_dist_para[0]

    dist_putative_para = dist_para + dist_alleles
    gps_overlap = []
    filter_seq = [(x[0], x[1]) for x in dist_putative_para if x[2] >= min_dist]

    for gp in gps:
        taxa_overlap = []
        for pair in itertools.combinations(gp, 2):
            item = tuple(
                sorted([pair[0][0] + "|" + pair[0][1], pair[1][0] + "|" + pair[1][1]])
            )
            if item in filter_seq:
                taxa_overlap.append(item)
        if taxa_overlap:
            gps_overlap.append(taxa_overlap)
    # Number of taxa that possess overlaping seq.
    n_taxa_overlap = len(
        list(set([x.split("|")[0] for x in list(itertools.chain(*filter_seq))]))
    )

    taxa_overlap_ratio = n_taxa_overlap / n_taxa_names

    return gps_overlap, taxa_overlap_ratio, dist_dict, max_dist


def check_paralogy_incompatibility_seq(tested_seq, component_content, paralog_pairs):
    """
    Function that scan the component to verify that adding
    a new sequence does not add a paralog to the component
    """
    result = True
    for item in component_content:
        if tuple(sorted([tested_seq, item])) in paralog_pairs:
            result = False
            break
    return result


def check_paralogy_incompatibility_comp(comp_1, comp_2, paralog_pairs):
    """ 
    Function that scan the component to verify that adding
    a new component does not add a paralog to the component
    """
    result = True
    possible_combinations = list(itertools.product(comp_1, comp_2))
    for item in possible_combinations:
        if tuple(sorted(item)) in paralog_pairs:
            result = False
            break
    return result


def obtain_dist_comps(
    fasta, taxa, min_overlap, taxa_profiles, gp_overlap, dist, dist_dict, MaxDist
):
    para_pairs = [
        tuple(sorted(x))
        for x in list(itertools.chain(*gp_overlap))
        if does_overlap(x, min_overlap, taxa_profiles)
    ]

    L_pairs = [
        [x[0], x[1], x[2]["distance"]]
        for x in dist
        if x[2]["distance"] < MaxDist
        and does_overlap((x[0], x[1]), min_overlap, taxa_profiles)
        and tuple(sorted([x[0], x[1]])) not in para_pairs
    ]

    sort_L_pairs = sorted(
        [x for x in L_pairs if x[0] in list(taxa.keys()) and x[1] in list(taxa.keys())],
        key=lambda x: x[2],
    )
    comp = {}
    #  taxa found in all components:
    taxa_all_comp = []
    for pair in sort_L_pairs:
        if pair[0] not in taxa_all_comp and pair[1] not in taxa_all_comp:
            taxa_all_comp.extend([pair[0], pair[1]])
            taxa_all_comp = list(set(taxa_all_comp))
            comp_id = max(list(comp.keys()) + [0]) + 1
            comp[comp_id] = [pair[0], pair[1]]
        elif pair[0] not in taxa_all_comp and pair[1] in taxa_all_comp:
            for k, v in comp.items():
                if pair[1] in v:
                    if check_paralogy_incompatibility_seq(pair[0], v, para_pairs):
                        comp[k].append(pair[0])
                        taxa_all_comp.append(pair[0])
                        taxa_all_comp = list(set(taxa_all_comp))
                        break
        elif pair[1] not in taxa_all_comp and pair[0] in taxa_all_comp:
            for k, v in comp.items():
                if pair[0] in v:
                    if check_paralogy_incompatibility_seq(pair[1], v, para_pairs):
                        comp[k].append(pair[1])
                        taxa_all_comp.append(pair[1])
                        taxa_all_comp = list(set(taxa_all_comp))
                        break
        else:
            pair0_key = None
            pair1_key = None
            for k, v in comp.items():
                if not pair0_key and pair[0] in v:
                    pair0_key = k
                if not pair1_key and pair[1] in v:
                    pair1_key = k
            if pair0_key != pair1_key:
                if check_paralogy_incompatibility_comp(
                    comp[pair0_key], comp[pair1_key], para_pairs
                ):
                    comp[pair0_key].extend(comp[pair1_key])
                    del comp[pair1_key]

    #  Remove components made of only two sequences that are identified as
    #  alleles (they contain the allele tag in their names)
    filter_comp = []
    for comp in comp.values():
        if len(comp) > 2:
            filter_comp.append(comp)
        elif comp[0].split("allele")[0] != comp[1].split("allele")[0]:
            filter_comp.append(comp)

    if filter_comp:
        return filter_comp
    else:
        return None


def get_para_score(pair, gps_overlap):
    """
    Function that yields the number of paralogs shared between two components
    for the common taxa
    """
    all_para_pairs = list(itertools.chain(*gps_overlap))
    shared_para = [
        sorted(x)
        for x in list(itertools.product(pair[0], pair[1]))
        if tuple(sorted(x)) in all_para_pairs
    ]
    taxa_para = set([x.split("|")[0] for x in list(itertools.chain(*shared_para))])
    #
    #   Out of the two components making up the pair,
    #   return the largest proportion calculated as the ratio of
    #   ((number of taxa involved into paralogy pairs) /
    #        (number of taxa contained in the component)
    #
    para_score = max(
        [
            len(taxa_para) / len(set([x.split("|")[0] for x in pair[0]])),
            len(taxa_para) / len(set([x.split("|")[0] for x in pair[1]])),
        ]
    )
    return para_score


def comp_incompatibility(tested_component, scrutinized_components, unlinked_components):
    """
    Function that check that linking two components does not create
    a cycle with a previously excluded component.
    Subfunction for cluster_comp_JI()
    """
    result = True
    for item in scrutinized_components:
        if (tested_component, item) in unlinked_components or (
            item,
            tested_component,
        ) in unlinked_components:
            result = False
            break
    return result


def set_comp_incompatibility(comp_1, comp_2, unlinked_components):
    """
    Function that check that linking two components does not create
    a cycle with a previously excluded component.
    Similar to comp_incompatibility() but this time one component and
    a set of components that is tested but two sets of components.
    Subfunction for cluster_comp_JI()
    """
    result = True
    possible_combinations = list(itertools.product(comp_1, comp_2))
    for item in possible_combinations:
        if item in unlinked_components or (item[1], item[0]) in unlinked_components:
            result = False
            break
    return result


def cluster_comp_JI(unlinked_comps, linked_comps):
    """
    Function that cluster components based on the smallest amout of paralogy
    """
    comps = {}
    cluster_comps = []
    for pair in linked_comps:
        if pair[0][0] not in cluster_comps and pair[0][1] not in cluster_comps:
            cluster_comps.extend([pair[0][0], pair[0][1]])
            cluster_comps = list(set(cluster_comps))
            comp_id = max(list(comps.keys()) + [0]) + 1
            comps[comp_id] = [pair[0][0], pair[0][1]]

        elif pair[0][0] not in cluster_comps and pair[0][1] in cluster_comps:
            for k, v in comps.items():
                if pair[0][1] in v:
                    if comp_incompatibility(pair[0][0], v, unlinked_comps):
                        comps[k].append(pair[0][0])
                        cluster_comps.append(pair[0][0])
                        cluster_comps = list(set(cluster_comps))
                        break

        elif pair[0][1] not in cluster_comps and pair[0][0] in cluster_comps:
            for k, v in comps.items():
                if pair[0][0] in v:
                    if comp_incompatibility(pair[0][1], v, unlinked_comps):
                        comps[k].append(pair[0][1])
                        cluster_comps.append(pair[0][1])
                        cluster_comps = list(set(cluster_comps))
                        break

        else:
            pair_0_0_key = None
            pair_0_1_key = None
            for k, v in comps.items():
                if not pair_0_0_key and pair[0][0] in v:
                    pair_0_0_key = k
                if not pair_0_1_key and pair[0][1] in v:
                    pair_0_1_key = k
            if pair_0_0_key != pair_0_1_key:
                if set_comp_incompatibility(
                    comps[pair_0_0_key], comps[pair_0_1_key], unlinked_comps
                ):
                    comps[pair_0_0_key].extend(comps[pair_0_1_key])
                    del comps[pair_0_1_key]
    return comps


def retrieve_comps_JI(jaccard_filter, gps_overlap, min_overlap_ratio):
    force_separation = True
    jaccard_para_score = {}
    for pair in jaccard_filter:
        para_score = get_para_score(pair[0], gps_overlap)
        jaccard_para_score[(tuple(pair[0][0]), tuple(pair[0][1]))] = para_score

    unlinked_comps = [
        k
        for k in list(jaccard_para_score.keys())
        if jaccard_para_score[k] >= min_overlap_ratio
    ]
    linked_comps = sorted(
        [
            (k, jaccard_para_score[k])
            for k in list(jaccard_para_score.keys())
            if jaccard_para_score[k] < min_overlap_ratio
        ],
        key=lambda x: x[1],
    )

    #  Cluster the components according to the amount of paralogy
    comps = cluster_comp_JI(unlinked_comps, linked_comps)

    clustered_items = set(list(itertools.chain(*list(comps.values()))))
    unclustered_items = set(
        list(itertools.chain(*list(jaccard_para_score.keys())))
    ).difference(clustered_items)
    #
    #  If clustering fails, the operation can be forced by using
    #  the largest of the paralogy distances as a threshold.
    #
    if len(comps) + len(unclustered_items) == 1 and force_separation is True:
        min_ratio_taxa_overlap = max(list(jaccard_para_score.values()))
        unlinked_comps = [
            k
            for k in list(jaccard_para_score.keys())
            if jaccard_para_score[k] >= min_ratio_taxa_overlap
        ]
        linked_comps = sorted(
            [
                (k, jaccard_para_score[k])
                for k in list(jaccard_para_score.keys())
                if jaccard_para_score[k] < min_ratio_taxa_overlap
            ],
            key=lambda x: x[1],
        )

        #  Cluster the components according to the amount of paralogy
        comps = cluster_comp_JI(unlinked_comps, linked_comps)
        clustered_items = set(list(itertools.chain(*list(comps.values()))))
        unclustered_items = set(
            list(itertools.chain(*list(jaccard_para_score.keys())))
        ).difference(clustered_items)

    if len(comps) + len(unclustered_items) == 1:
        final_comps = []
        final_comps.append(list(itertools.chain(*[list(x) for x in comps.values()[0]])))
        return final_comps
    else:
        all_clustered_comps = []
        for v in comps.values():
            all_clustered_comps.append(list(itertools.chain(*[list(x) for x in v])))
        return all_clustered_comps + [list(x) for x in unclustered_items]


def inparalog_separation(fasta, cfg, **kargs):
    """
    Main function for STEP5:
    interp. 
        seq_overlap_splitting (int)     The minimum overlap between two paralogs used
                                        to compute the minimum distance between paralogs
                                        that allows to split the distance graph into
                                        connected components.
        seq_overlap_merging (int)       The minimum overlap between two sequences  
                                        required for using pairwise distances for 
                                        adding sequences to components.
        outgroups (list(str))           List of the outgroups used to separate putative
                                        outparalogs.
        max_jaccard_value (float)       Value used during the component fusion step.
                                        Pairs of components that exceed this value of 
                                        theJaccard Index for their taxa are not 
                                        considered for fusion. Lower value means higher
                                        stringency.
        max_alleles_distance (float)    indicates the max distance allowed between 2 
                                        alleles.
        paralog_allele_distance_ratio (float) ratio between minimum distance between
                                        paralogues and maximum distance between alleles
        large_component_ratio (float)   During the component clustering step, large
                                        components possess a number of taxa larger than 
                                        'large_component_ratio'*#taxa largest component.
                                        These components are tested for union, smaller
                                        components are desassembled and their sequences
                                        are assigned with distances to the largest 
                                        components.
        min_taxa_in_alignment (int)     Number of taxa to serve as threshold between
                                        what is deemed to be a 'complete' alignment and  
                                        an 'imcomplete' one.                                                                      
    """
    in_suffix = cfg["input_suffix"]
    in_format = cfg["input_format"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    min_overlap = cfg["seq_overlap_merging"]
    min_overlap_ratio = cfg["min_proportion_taxa_that_overlap"]
    excluded_para = cfg["excluded_paralogs"]
    # max_dist = cfg["max_distance"]
    dist_ratio = cfg["paralog_allele_distance_ratio"]
    large_comp_ratio = cfg["large_component_ratio"]
    outgroups = cfg["outgroups"]
    max_jaccard = cfg["max_jaccard_value"]
    min_taxa = cfg["min_taxa_in_alignment"]
    overlap_max = cfg["seq_overlap_splitting"]
    max_alleles_dist = cfg["max_alleles_distance"]
    basename = fasta.split(in_format)[0]

    error_template = "An exception of type {0} occurred when trying to \
 {2}. Arguments:\n{1!r}"
    DEBUG = True

    try:
        assert (
            0 <= cfg["max_jaccard_value"] <= 1
            and 0 <= cfg["large_component_ratio"] <= 1
        )
    except:
        s = "Wrong parameters: 'max_jaccard_value' \
         and large_component_ratio must be in [0, 1]"
        raise ValueError(s)
        sys.exit(s)

    try:
        cfg["dna_model"] in [0, 1, 2, 3, 4, 5]
    except:
        s = "Wrong parameters: 'dna_model' must be in [0, 1, 2, 3, 4, 5]"
        raise ValueError(s)
        sys.exit(s)

    try:
        records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
    except:
        exception = "[Error]: Cannot open alignment {}\n\
        Verify that the file is in fasta format".format(
            fasta
        )
        return exception

    n_seq = len(records)
    if n_seq < min_taxa:
        s_i = "STEP5 Done with fasta: {}".format(fasta)
        s_d = "No transformation performed on fasta {0}\n\
        which only contains {1} sequences".format(
            fasta, n_seq
        )
        return (s_i, s_d)

    length = len(list(records.values())[0])
    try:
        assert min_overlap < length and min_overlap_ratio < length
    except:
        exception = "[Error]: for fasta {} 'seq_overlap_splitting'\
        or 'seq_overlap_merging' exceeds the length of the alignment".format(
            fasta
        )
        return (exception, None)

    #  1- Remove empty sequences
    for record in records.keys():
        if set(list(str(records[record].seq))) == set(["-"]):
            records.pop(record, None)

    #  Calculate distance matrix used to assign seq to comp
    calculate_distance_matrix(
        fasta, {"input_suffix": in_format, "dna_model": cfg["dna_model"]}
    )
    taxa = {records[rec].name: str(records[rec].seq) for rec in records.keys()}

    "for each sequence we recode 0 for '-', 1 otherwise,\
    into a dictionary {name_number: profile}"
    taxa_profiles = sequence_profiler(taxa)

    distance_matrix = basename + "_for_dnadist.fasta.distmat"
    dist = distance_matrix_parser(distance_matrix)

    taxa_no_outgp_para = {
        x: y
        for x, y in taxa.items()
        if x.split("|")[0] not in outgroups and x.split("|")[0] not in excluded_para
    }

    taxa_no_outgp = {x: y for x, y in taxa.items() if x.split("|")[0] not in outgroups}

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
                get_overlap_score(tuple(sorted([x, seq])), taxa_profiles) > overlap_max
            )
        ]
    inpara_validation_input = [
        taxa_no_outgp,
        taxa_profiles,
        min_overlap,
        taxa_no_outgp_para,
        dist,
        dist_ratio,
        dict_overlap_seq,
        max_alleles_dist,
    ]

    #  Retrieve the empirical max distance between paralogs: MaxDist
    inpara_validation_output = inpara_validation(*inpara_validation_input)
    gp_overlap, taxa_overlap, dist_dict, MaxDist = inpara_validation_output

    # 1- No max dist can be calculating for separating paralogs.
    #    Return the original alignment.
    if MaxDist is None:
        SeqIO.write(
            [x[1] for x in sorted(records.items(), key=lambda x: x[0])],
            basename + out_suffix + "_no_threshold.fasta",
            "fasta",
        )
        shutil.move(basename + out_suffix + "_no_threshold.fasta", out_folder)
        s_i = "Done with fasta {}".format(fasta)
        s_d = "No distance threshold between paralogs can be \
calculated for {} Returning original alignment".format(
            fasta
        )
        return (s_i, s_d)

    #  2- Max dist is calculated but sequences fail to sufficiently overlap.
    #     Paralogy is not detected. Return the original alignment.
    elif taxa_overlap < min_overlap_ratio:
        print(
            "Sequence overlap below the set threshold. \n\
        {} \n\is considered inparalog free".format(
                fasta
            )
        )
        SeqIO.write(
            [x[1] for x in sorted(records.items(), key=lambda x: x[0])],
            basename + out_suffix + "_no_overlap.fasta",
            "fasta",
        )
        shutil.move(basename + out_suffix + "_no_overlap.fasta", out_folder)
        s_i = "Done with fasta {}".format(fasta)
        s_d = "Sequence overlap below the set threshold. \n\
{} is considered inparalog free".format(
            fasta
        )
        return (s_i, s_d)

    #  3- Max dist is calculated and sequences do sufficiently overlap.
    #     Paralogy is detected. Attempt a clustering based on dist method
    else:
        #  Switch that records a successful separation of components
        separation = False
        if DEBUG:
            print("STEP 3: Attempting a  distance based clustering")
        # Compute distance based components"
        comp_1 = obtain_dist_comps(
            fasta,
            taxa_no_outgp,
            min_overlap,
            taxa_profiles,
            gp_overlap,
            dist,
            dist_dict,
            MaxDist,
        )

        #  Find large components
        if comp_1:
            max_comp = max([len(x) for x in comp_1])
            large_comps = [c for c in comp_1 if len(c) > (large_comp_ratio * max_comp)]

        #  3.1- The component list is not empty and contains at least
        #       two large components.
        if comp_1 is not None and len(large_comps) > 1:
            jaccard = []
            for pair in itertools.combinations(large_comps, 2):
                jaccard.append(calculate_Jaccard_index(pair))
            jaccard_filtered = [x for x in jaccard if x[1] <= max_jaccard]

            #  case 1: Components cannot be clustered by Jaccard.
            #          Separation is  succesful.
            if len(jaccard_filtered) < 2:
                comp_2 = large_comps
                separation = True
            #  case 2: Components are sufficiently different to test
            #           for components clustering
            else:
                comp_2 = retrieve_comps_JI(
                    jaccard_filtered, gp_overlap, min_overlap_ratio
                )
                if len(comp_2) > 1:
                    separation = True

        #  3.2- The separation attempt into several components was
        #       unsuccessful. Return the original alignment.
        if not separation:
            SeqIO.write(
                [x[1] for x in sorted(records.items(), key=lambda x: x[0])],
                basename + out_suffix + "_no_separation.fasta",
                "fasta",
            )
            shutil.move(basename + out_suffix + "_no_separation.fasta", out_folder)
            s_i = "Done with fasta {}".format(fasta)
            s_d = "No component separation was obtained for fasta\n\ {}\n\
            Returning original alignment".format(
                fasta
            )
            return (s_i, s_d)

        #    3.3- Some components have been isolated
        #       Adding the remaining sequences using distances
        elif separation:
            #  Common parameter for all calls of dist_assignment()
            dist_assignmet_para = [
                basename,
                out_folder,
                MaxDist,
                min_taxa,
                outgroups,
                min_overlap,
                out_suffix,
            ]

            final_comps = dist_assignment(
                comp_2, records, taxa_profiles, dist, *dist_assignmet_para
            )

            write_fasta(
                final_comps, records, basename, out_folder, min_taxa, out_suffix
            )
            s_i = "Done with fasta {}".format(fasta)
            s_d = "Component separation was succesful for fasta {0}\n\
            Returning {1} alignements".format(
                fasta, len(final_comps)
            )
            return (s_i, s_d)



