#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Tools that separates outparalogs.
    1. Alleles that have been phased are identified by their 'allele' tag name.
        The longest allele is kept for further analyses.
    2. Remove empty sequences.
    3. Return the original alignment if none of the specified outgroups is present.
    4. Kmer pairwise identity is used to compute a kmer score. 
        For a pair of sequence, the Kmer score is calculated as:
        #positions with identity/total #positions where there is sequence overlap. 
        Different Kmer are tested in order to achieve the best separation.
    5. A weighted undirected graph is built with sequences as vertices and Kmer scores 
        as edges.
        The maximum Kmer score (max score) between an ingroup and an outgroup sequence 
        computed. All edges longer than max score are removed, yielding connected 
        commponents. If a single component is obtained, return the original alignment.
    6. When several components have been obtained, we test for overseparation by 
        computing the the Jaccard index for the number of taxa in the components:
        J(C1,C2) = (taxa in C1 and C2)/(taxa in C1 or C2). If the index's value
        is below a present threshold the sequence data is check to find if it supports
        the putative association based on the 'min_association_ratio' parameter.
        E.g. for a parameter value of 0,75 and with C1 that contains sequences (S0, S1) 
        and C2 (S2, S3), if based on Kmer identity the ratio between the number of 
        positions that support the grouping of C1 and C2 and the number of positions
        where there is sequence overlap between the two componentes is higher than 0.75.
        Also a component cannot have more than 'max_allele' overlapping sequences from
        the same individual.
    7. Distances (corrected or uncorrected according to user selection) are computed 
        between each sequence pair using the EMBOSS distmat software.
    8. If there is more than two components or two components and ungroup sequences,
        the smallest components (less than half the size of the largest) are 
        desassembled. Their sequences and the ungroup sequences are assigned to the 
        larger components based on smaller distance to one of the sequences in the 
        components distance as long as the distance are below the maximum distance
        threshold set by the user.
    9. Files that harbour less taxa than the set threshold are tagged as '_incomplete',
        otherwise they are tagged as '_complete'. Sequences that remained unassigned
        after the distance clustering are labelled as 'unassigned_seq'.
"""
__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import shutil
import operator
from collections import Counter
import networkx as nx
import itertools

from Bio import SeqIO

try:
    from aligera_scripts.utilities import (
        calculate_distance_matrix,
        distance_matrix_parser,
        sequence_profiler,
        remove_ambiguities,
        run_subprocess,
    )
except:
    from utilities import (
        calculate_distance_matrix,
        distance_matrix_parser,
        sequence_profiler,
        remove_ambiguities,
        run_subprocess,
    )


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def validate_associations(jaccard_filtered_pair, aln_dict_values):
    """
    Check whether the possible associations proposed by the jaccard index
    are actually found in the alignment.
    """

    #  jaccard_filtered_pair contains two networkx.classes.reportviews.NodeView objects
    #  that need to be converted to lists in order to be added to each other.
    pair = tuple(list(jaccard_filtered_pair[0][0]) + list(jaccard_filtered_pair[0][1]))
    #  positions where the association is not rejected,
    #  i.e. the sequences are not part of several components
    positive_positions = 0
    #  positions that reject the association,
    #  i.e. the sequences are divided between several components
    negative_positions = 0
    for position in aln_dict_values:
        group_intersections = sum(
            [1 if set(pair).intersection(set(x)) != set() else 0 for x in position]
        )
        if group_intersections == 1:
            positive_positions += 1
        elif group_intersections > 1:
            negative_positions += 1
    return (positive_positions, negative_positions)


def write_fasta(final_components, records, basename, out_folder, min_taxa, out_suffix):
    """
    Write the sequences in the final components to fasta.
    A complete align contains more sequences than the defined minimum,
     an incomplete otherwise. Final_components is a dict
     {(name_1,name_2): [name_1, name_2]}
    """
    i = 0
    j = 0
    for comp in final_components.values():
        if len(list(set([x.split("|")[0] for x in comp]))) >= min_taxa:
            name = basename + out_suffix + "_complete_{}.fasta".format(i)
            SeqIO.write([records[x] for x in sorted(comp)], name, "fasta")
            try:
                shutil.move(name, out_folder)
            except:
                exception = "[Error]: Cannot move file {0} to folder {1}".format(
                    name, out_folder
                )
                return exception
            i += 1
        else:
            name = basename + out_suffix + "_incomplete_{}.fasta".format(j)
            SeqIO.write([records[x] for x in sorted(comp)], name, "fasta")
            try:
                shutil.move(name, out_folder)
            except:
                exception = "[Error]: Cannot move file {0} to folder {1}".format(
                    name, out_folder
                )
                return exception
            j += 1


def dist_assignment(
    components,
    records,
    taxa_profiles,
    basename,
    out_folder,
    max_dist,
    min_taxa,
    outgroups,
    min_overlap,
    out_suffix,
    large_comp_ratio,
):
    """
    Large and small components are separated.
    Sequences in small components are reassigned based on distance
    Threshold between small and larged is set to 
    (size large component) * large_comp_ratio.
    """
    comps = components
    max_length_component = sorted(comps, key=lambda x: len(x))[-1]
    large_comps = [
        c for c in comps if len(c) > len(max_length_component) * large_comp_ratio
    ]
    final_comps = {tuple(comp): comp for comp in large_comps}
    unassigned_seqs = list(
        set(list(records.keys())).difference(set(itertools.chain(*large_comps)))
    )
    #
    #  If there is at least one large components and some small components
    #
    if len(large_comps) >= 1 and len(unassigned_seqs) > 0:
        distance_matrix = basename + "_for_dnadist.fasta.distmat"
        distances = distance_matrix_parser(distance_matrix)
        if not distances:
            exception = "[Error] problem parsing distances\nExiting"
            raise Exception(exception)
            sys.exit()
        unassignable_taxa = []
        for unassigned_seq in unassigned_seqs:
            min_dist = {}
            for comp in large_comps:
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
            #
            #  If one of the components already contains the unassigned taxon,
            #  only a limited overlap between the two sequences is allowed
            #
            if sorted_min_dist:
                #
                #  List of pairs containing one unassigned seq and
                #  one seq from the selected component
                #
                overlaps = [[unassigned_seq, x] for x in sorted_min_dist[0][0]]
                #
                #  Calculate overlap between unassigned seq and each member
                #  of component. Return the max.
                #
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
                #
                #  Case 1. There is overlap with seq from the same taxon and
                #   the dist is smaller than the threshold
                #   or the overlap is small but the distance is minute
                #

                if unassigned_seq.split("|")[0] not in outgroups:
                    if (
                        sorted_min_dist[0][1] <= max_dist and max_overlap > min_overlap
                    ) or sorted_min_dist[0][1] <= max_dist / 3:
                        new_final_comps = list(final_comps[sorted_min_dist[0][0]])
                        new_final_comps.append(unassigned_seq)
                        final_comps[sorted_min_dist[0][0]] = new_final_comps
                    else:
                        unassignable_taxa.append(unassigned_seq)
                #
                #  Case 2. We are dealing with the outgroup and
                #   we double the max dist for assigning seq to component
                #
                else:
                    if (
                        sorted_min_dist[0][1] <= max_dist * 2
                        and max_overlap > min_overlap
                    ) or sorted_min_dist[0][1] <= max_dist / 3:
                        new_final_comps = list(final_comps[sorted_min_dist[0][0]])
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


def calculate_Jaccard_index(pair):
    """
    Compute the Jaccard index from taxa present in a pair of components
    pair is a tupple(component1, component2)
    """
    species_1 = [x.split("|")[0] for x in pair[0]]
    species_2 = [x.split("|")[0] for x in pair[1]]
    JI = float(len(set(species_2).intersection(set(species_1)))) / len(
        set(species_1).union(set(species_2))
    )
    return (pair, JI)


def get_alleles(taxa, outgroups):
    """
    Identify the alleles and remove the shortest allele
    The alleles are identified if they correspond to
    the allele_0/allele_1 designation obtained from the phasing step.
    taxa is a {name: sequence} dict.
    """
    taxa_with_alleles = list(
        set(
            [
                x[0]
                for x in [
                    y.split("allele")
                    for y in taxa.keys()
                    if y.split("|")[0] not in outgroups
                ]
                if x[0] + "allele_0" in list(taxa.keys())
                and x[0] + "allele_1" in list(taxa.keys())
            ]
        )
    )
    alleles = [(p + "allele_0", p + "allele_1") for p in taxa_with_alleles]

    sorted_alleles = sorted(
        [
            (
                (allele_0, len(remove_ambiguities(taxa[allele_0], "Nn-"))),
                (allele_1, len(remove_ambiguities(taxa[allele_1], "Nn-"))),
            )
            for allele_0, allele_1 in alleles
        ],
        key=lambda x: x[1],
    )
    longest_alleles = [x[1][0] for x in sorted_alleles]
    return longest_alleles


def get_target_lists(taxa, aln_length, kmer, outgroups):
    """
    For a given kmer several lists are returned:
    seq_group: dict(idx: components list), {2:[[1,5],[2,4]],3:[[1,5],[2,3,4]]}
    pairs_to_compare_with_outgroup: list pairs containing a single outgroup seq
    pairs_to_compare_without_outgroup: list pairs without outgroup seq
    list_names_per_position: at each position return a list of sequence names:
        [[1,2,4,5],[1,2,3,4,5]]
    """
    #   Identify the alleles and remove the shortest allele.
    longest_alleles = get_alleles(taxa, outgroups)
    #  Names of all sequences that do not have alleles
    #  and names with the sequences that do have an allele pair.
    seqs_group = [x for x in taxa.keys() if "allele" not in x]
    seqs_group.extend(longest_alleles)
    #
    #  A dictionary of dictionaries, keys are position values are subdict.
    #  In the subdict, keys are sequence names and values are kmer
    #  aln_dict.items() >>> [(2, {1: 'TTT'}), (5, {1: 'AAT'}), (6, {0: 'CCA', 1: 'ATT'})
    aln_dict = {}
    for name in sorted(seqs_group):
        s_kmers = list(
            enumerate(
                [
                    tuple([taxa[name][idx : idx + kmer], name])
                    for idx in range(aln_length)
                ]
            )
        )
        for item in s_kmers:
            if set(item[1][0]).intersection(set(["-", "N", "n"])) == set():
                if item[0] in list(aln_dict.keys()):
                    aln_dict[item[0]][item[1][1]] = item[1][0]
                else:
                    aln_dict[item[0]] = {item[1][1]: item[1][0]}

    def keys_with_same_value_as_key(d, k):
        return tuple([key for key in d.keys() if d[key] == d[k]])

    #  seq_group_dict.items() >>>[(2, [[1,5],[2,3,4]]), (3, [[1,2,5],[2,3,4]])]
    #  with keys as posisions and sequences names grouped into components
    #  larger than 1
    seq_group_dict = {}
    for idx in aln_dict.keys():
        d = aln_dict[idx]
        L = [
            keys_with_same_value_as_key(d, k)
            for k in d.keys()
            if len(keys_with_same_value_as_key(d, k)) > 1
        ]
        seq_group_dict[idx] = [list(x) for x in set(L)]

    #  Counter that yields the number of component each sequence is part of
    flat_seq_group_dict_values = []
    for idx in seq_group_dict.values():
        for component in idx:
            flat_seq_group_dict_values.extend(component)
    components_per_name_counter = Counter(flat_seq_group_dict_values)
    #
    #  Check all pairs present in components.
    #  Each component is recoded as a clique
    list_pairs_in_components = set()
    #  for each position all the pairs are recorded
    list_pairs_per_position = []
    #  sequence names present at each position
    list_names_per_position = []

    for idx in seq_group_dict.values():
        clique_per_position = []
        names_per_position = []
        for component in idx:
            names_per_position.extend(component)
            clique_per_position.append(component)
            cliques = [
                tuple(sorted(x)) for x in list(itertools.combinations(component, 2))
            ]
            list_pairs_in_components = list_pairs_in_components.union(set(cliques))

        list_pairs_per_position.append(clique_per_position)
        list_names_per_position.append(names_per_position)

    #  A set of all pairs present among all components.
    #  These are the pairs that will be used for distance calculations
    set_pairs_in_components = list(set(list_pairs_in_components))

    #  pairs_to_compare contains taxa that are present in many components
    pairs_to_compare = [
        x
        for x in set_pairs_in_components
        if components_per_name_counter[x[0]] > 100
        and components_per_name_counter[x[1]] > 100
    ]

    outgroup_seqs = [x for x in taxa.keys() if x.split("|")[0] in outgroups]

    pairs_to_compare_with_outgroup = [
        x
        for x in pairs_to_compare
        if (
            (x[0] in outgroup_seqs and x[1] not in outgroup_seqs)
            or (x[1] in outgroup_seqs and x[0] not in outgroup_seqs)
        )
    ]

    pairs_to_compare_without_outgroup = [
        x
        for x in pairs_to_compare
        if (
            x not in pairs_to_compare_with_outgroup
            and x[0] not in outgroup_seqs
            and x[1] not in outgroup_seqs
        )
    ]

    return (
        seq_group_dict,
        pairs_to_compare_with_outgroup,
        pairs_to_compare_without_outgroup,
        list_pairs_per_position,
        list_names_per_position,
    )


def compute_component_score(
    pair, list_pairs_per_position, list_names_per_position, min_overlap=100
):
    """
    together + 1 when a sequence pair is in the same component at a given position 
    using Kmer identity.
    separated + 1 otherwise.       
    """
    together = 0
    separated = 0
    for idx, components in enumerate(list_pairs_per_position):
        if (
            pair[0] in list_names_per_position[idx]
            and pair[1] in list_names_per_position[idx]
        ):
            switch = False
            for component in components:
                if pair[0] in component and pair[1] in component:
                    together += 1
                    switch = True
                    break
            if not switch:
                separated += 1
    if together and (together + separated) > min_overlap:
        return (pair, float(together) / (together + separated), together + separated)
    else:
        return (pair, 0.0, together + separated)


def does_overlap(item, min_overlap, taxa_profiles):
    """
    Calculte sequence overlap
    """
    item_overlap = sum(
        a == b for a, b in zip(taxa_profiles[item[0]], taxa_profiles[item[1]]) if a == 1
    )
    if item_overlap > min_overlap:
        return True
    else:
        return False


def obtain_components(
    fasta, taxa, idxs_length, kmer, outgroups, min_overlap, taxa_profiles
):
    target_lists = get_target_lists(taxa, idxs_length, kmer, outgroups)
    #  dictionary: key are position, values are components at the position
    seq_group_dict = target_lists[0]
    #  list of empirical pairs that have one of the sequences from the outgroup
    pairs_with_outgroup = target_lists[1]
    #  list of empirical pairs that have no outgroup sequence
    pairs_without_outgroup = target_lists[2]
    #  list of empirtical pairs at each position
    pairs_per_idx = target_lists[3]
    #  list of sequence names present at each position
    names_per_idx = target_lists[4]
    #  Filter the empirical pair according to their overlap
    pairs_outgroup_filtered = [
        x for x in pairs_with_outgroup if does_overlap(x, min_overlap, taxa_profiles)
    ]
    pairs_no_outgroup_filtered = [
        x for x in pairs_without_outgroup if does_overlap(x, min_overlap, taxa_profiles)
    ]
    #
    #  Try to get the best score for the maximum difference with outgroup
    #  by testing several values for min_overlap
    dist_ougroup = [
        compute_component_score(pair, pairs_per_idx, names_per_idx)
        for pair in pairs_outgroup_filtered
    ]

    #  Attempt to calculate the maximum distance to the outgroup
    #  Several min overlap are tested from highest to lowest as
    #   a high overlap is more precise but harder to obtain than
    #   a low overlap
    for min_seq_overlap in [300, 200, 100]:
        sorted_dist_ougroup = sorted(
            [(x[0], x[1]) for x in dist_ougroup if x[2] > min_seq_overlap],
            key=lambda x: x[1],
        )
        if sorted_dist_ougroup:
            max_dist_ougroup = sorted_dist_ougroup[-1][1]
            selected_outgroup = [
                x.split("|")[0]
                for x in sorted_dist_ougroup[-1][0]
                if x.split("|")[0] in outgroups
            ][0]
            break
        else:
            max_dist_ougroup = None

    if not max_dist_ougroup:
        return (None, None, None)

    dist_no_outgroup = [
        compute_component_score(pair, pairs_per_idx, names_per_idx)
        for pair in pairs_no_outgroup_filtered
    ]
    if not dist_no_outgroup:
        return (None, None, None)

    G = nx.Graph()
    
    #  Filter away distances that are larger than the threshold set by
    #   the maximum distance to the outgroup and are calculated from
    #   sufficient sequence overlap
    
    G.add_edges_from(
        [
            item[0]
            for item in dist_no_outgroup
            if item[1] > max_dist_ougroup and item[2] > min_overlap
        ]
    )

    components = [
        component.nodes() for component in nx.connected_component_subgraphs(G)
    ]

    return components, seq_group_dict, selected_outgroup


def validate_component(comp, taxa_profiles, min_overlap, max_alleles):
    """
    A component is validated if for the same taxon it does not contain
    more overlapping alleles than the set maximum.
    """
    group_by_taxa = []
    for key, gp in itertools.groupby(
        sorted([x.split("|") for x in comp], key=lambda x: x[0]), lambda x: x[0]
    ):
        group_by_taxa.append(list(gp))

    groups = [x for x in group_by_taxa if len(x) > 1]
    count_taxa_with_overlap = 0
    for gp in groups:
        for pair in itertools.combinations(gp, 2):
            x = tuple([pair[0][0] + "|" + pair[0][1], pair[1][0] + "|" + pair[1][1]])
            if does_overlap(x, min_overlap, taxa_profiles):
                count_taxa_with_overlap += 1
                break
    if count_taxa_with_overlap < max_alleles:
        return True
    return False


def outparalog_separation(fasta, cfg, **kargs):
    """
    Main function for STEP4:
    interp. 
        kmer (int)                      The number of consecutive bases to be examine
                                        for computing the Kmer score.
        min_sequences_overlap (int)     The minimum overlap between two sequences  
                                        required for computing pairwise distances with
                                        the EMBOSS dismat software.
        outgroups (list(str))           List of the outgroups used to separate putative
                                        outparalogs.
        max_alleles (int)               Maximum number of putative alleles allowed 
                                        within a cluster during the cluster fusion step. 
                                        The final number of putative alleles in the 
                                        alignemnt can exceed this value in latter steps.
        max_jaccard_value (float)       Value used during the component fusion step.
                                        Pairs of components that exceed this value of 
                                        the Jaccard Index (JI) for their taxa are not 
                                        considered for fusion. Lower value means higher
                                        stringency.
        min_association_ratio (float)   Parameter that together with the Jaccard Index
                                        value controles the stringency of the fusion of
                                        components. A higher ratio increases the 
                                        stringency.
        max_distance (float)            Max distance between sequences allowed during
                                        the assignment of unclustered sequences to their
                                        components.
        min_taxa_in_alignment (int)     Number of taxa to serve as threshold between
                                        what is deemed to be a 'complete' alignment and  
                                        an 'imcomplete' one.
        large_component_ratio (float)   Value of the Size ratio of the largest component
                                        that is used to distinguish small components:
                                        a small component has
                                        size < (size largest) * large_component_ratio.
                                        Large components are used during the clustering
                                        phase with JI.
    """
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    kmers = cfg["kmers"]
    min_overlap = cfg["min_sequences_overlap"]
    outgroups = cfg["outgroups"]
    max_alleles = cfg["max_alleles"]
    max_jaccard = cfg["max_jaccard_value"]
    association_ratio = cfg["min_association_ratio"]
    max_dist = cfg["max_distance"]
    min_taxa = cfg["min_taxa_in_alignment"]
    large_comp_ratio = cfg["large_component_ratio"]
    basename = fasta.split(in_suffix)[0]
    error_template = "An exception of type {0} occurred when trying to \
 {2}. Arguments:\n{1!r}"
    DEBUG = True

    try:
        assert (
            0 <= cfg["min_association_ratio"] <= 1
            and 0 <= cfg["max_jaccard_value"] <= 1
            and 0 <= cfg["large_component_ratio"] <= 1
        )
    except:
        s = "Wrong parameters: 'min_association_ratio', 'max_jaccard_value' \
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
    try:
        assert n_seq > 3
    except:
        s_i = "STEP4 Done with fasta: {}".format(fasta)
        s_d = "No transformation performed on fasta {0}\n\
        which only contains {1} sequences".format(
            fasta, n_seq
        )
        return (s_i, s_d)

    length = len(list(records.values())[0])
    try:
        assert min_overlap < length
    except:
        exception = "[Error]: for fasta {} 'min_sequences_overlap'\
        exceeds the length of the alignment".format(
            fasta
        )
        return (exception, None)
    if DEBUG:
        print("working on fasta {}".format(fasta))
    #  1- Remove empty sequences.
    if DEBUG:
        print("1- Remove empty sequences.")
    for record in records.keys():
        if set(list(str(records[record].seq))) == set(["-"]):
            records.pop(record, None)

    #  2- Return the original aln if no outgroup
    if (
        set(outgroups).intersection(set([x.split("|")[0] for x in records.keys()]))
        == set()
    ):
        SeqIO.write(
            [x[1] for x in sorted(records.items(), key=lambda x: x[0])],
            basename + out_suffix + "_no_outgroup.fasta",
            "fasta",
        )
        shutil.move(basename + out_suffix + "_no_outgroup.fasta", out_folder)
        s_i = "Done with fasta {}".format(fasta)
        s_d = "No outgroup is present for fasta {}".format(fasta)
        return (s_i, s_d)

    #  Calculate distance matrix used to assign sequences to components
    if DEBUG:
        print("Calculate distance matrix")
    calculate_distance_matrix(fasta, cfg)
    taxa = {records[rec].name: str(records[rec].seq) for rec in records.keys()}
    #  For each sequence we recode 0 for '-', 1 otherwise,
    #  into a dictionary {name_number: profile}
    taxa_profiles = sequence_profiler(taxa)
    if DEBUG:
        print("Test different kmer in the kmers list")
    #  3- Test different kmer in the kmers list in order to obtain several components.
    #
    retained_comp = None
    for kmer in kmers:
        idxs_length = len(list(records.values())[0].seq) - kmer
        if idxs_length <= 20:
            s_i = "Done with fasta {}".format(fasta)
            s_d = "In fasta {0} with kmer {1}, \n\
            the alignment is too short to perform STEP4".format(
                fasta, kmer
            )
            return (s_i, s_d)

        comp_1, aln_dict, accepted_outgroup = obtain_components(
            fasta, taxa, idxs_length, kmer, outgroups, min_overlap, taxa_profiles
        )
        if DEBUG:
            if comp_1:
                print("# components: {0} with Kmer {1}".format(len(comp_1), kmer))
        if comp_1 is None:
            continue
        elif comp_1 is not None and len(comp_1) == 1:
            if validate_component(comp_1[0], taxa_profiles, min_overlap, max_alleles):
                retained_comp = comp_1
                break
            else:
                retained_comp = comp_1
        else:
            retained_comp = comp_1
            break

    #  4- If no components are obtained, return the original aln.
    #
    if retained_comp is None or len(retained_comp) == 0:
        SeqIO.write(
            [x[1] for x in sorted(records.items(), key=lambda x: x[0])],
            basename + out_suffix + "_no_splitting.fasta",
            "fasta",
        )
        try:
            shutil.move(basename + out_suffix + "_no_splitting.fasta", out_folder)
        except Exception as ex:
            s = "move fasta: {}".format(fasta)
            message = error_template.format(type(ex).__name__, ex.args, s)
            print(message)

        s_i = "Done with fasta {}".format(fasta)
        s_d = "No components were found,\n\
            using kmer {}: no inparalog removal was attempted".format(
            kmer
        )
        return (s_i, s_d)

    #  4. If a single component is found, use it to assign the sequence
    #     no included based on a distance method.
    if DEBUG:
        print("Cluster with distances")
    #  Common parameter for all calls of dist_assignment()
    dist_assignmet_para = [
        basename,
        out_folder,
        max_dist,
        min_taxa,
        outgroups,
        min_overlap,
        out_suffix,
        large_comp_ratio,
    ]
    if len(retained_comp) == 1:
        final_comp = dist_assignment(
            retained_comp, records, taxa_profiles, *dist_assignmet_para
        )
        if DEBUG:
            print("writing a single component")
        write_fasta(final_comp, records, basename, out_folder, min_taxa, out_suffix)

        s_i = "Done with fasta {}".format(fasta)
        s_d = "# retained components: {}".format(len(final_comp))
        return (s_i, s_d)

    #  5. If there are some components, compute the Jaccard index
    #            between pair of components:
    else:
        jaccard = []
        for pair in itertools.combinations(retained_comp, 2):
            jaccard.append(calculate_Jaccard_index(pair))
        jaccard_filtered = [x for x in jaccard if x[1] <= max_jaccard] ###
        #
        #  6. case 1: Jaccard index can be calculated,
        #     but is set too high therefore clustering all components.
        #     No additional component clustering is performed and seq
        #     are assigned based on a distance method.
        if len(jaccard_filtered) == 0:
            if DEBUG:
                print("Cluster without Jaccard")
            final_comp = dist_assignment(
                comp_1, records, taxa_profiles, *dist_assignmet_para
            )
            if DEBUG:
                print("writing a several component")
            write_fasta(final_comp, records, basename, out_folder, min_taxa, out_suffix)

            s_i = "Done with fasta {}".format(fasta)
            s_d = "# retained components: {}".format(len(comp_1))
            return (s_i, s_d)
        #
        #  6. case 2: Jaccard index can be calculated and is below the set
        #     maximum. Associations between components are checked to see
        #     if they can be justified given the data.
        #     Remaining individual sequences are clustered by distance.
        elif len(jaccard_filtered) > 0:
            #  At each position we calculate the components.
            if DEBUG:
                print("Cluster with Jaccard")
            JI_gp = {}
            for pair in jaccard_filtered:
                gp = validate_associations(pair, list(aln_dict.values()))
                JI_gp[(tuple(pair[0][0]), tuple(pair[0][1]))] = {
                    "JI": pair[1],
                    "gp": gp,
                }

            edges = [
                (
                    key[0],
                    key[1],
                    {
                        "dist": JI_gp[key]["gp"][0]
                        / (JI_gp[key]["gp"][0] + JI_gp[key]["gp"][1])
                    },
                )
                if (JI_gp[key]["gp"][1] + JI_gp[key]["gp"][0]) != 0
                else (key[0], key[1], {"dist": 0.0})
                for key in JI_gp.keys()
            ]
            filtered_edges = [e for e in edges if e[2]["dist"] > association_ratio]

            #  Separate components clustered with Jaccard index
            #   from unclustered components.
            single_edges = []
            for e in filtered_edges:
                single_edges.append(e[0])
                single_edges.append(e[1])
            #  Components that are linked by an edge
            L_single_edges = list(set(single_edges))
            #  Components that lack an edge to any other component
            #  are separated from the rest
            unclustered_comp = []
            for e in edges:
                if e[0] not in L_single_edges:
                    unclustered_comp.append(e[0])
                if e[1] not in L_single_edges:
                    unclustered_comp.append(e[1])
            L_unclustered_comp = list(set(unclustered_comp))
            #  Compatible components are clustered and the remaining
            #  components are assigned using a distance method.

            if filtered_edges:
                #  Assembly of components that are associated with JI
                G_2 = nx.Graph()
                G_2.add_edges_from(filtered_edges)
                filtered_comp = [
                    list(set(itertools.chain(*c.nodes())))
                    for c in nx.connected_component_subgraphs(G_2)
                ]
                comp_2 = filtered_comp + [list(x) for x in L_unclustered_comp]
                #  add the seq from components that have not been clustered
                final_comp = dist_assignment(
                    comp_2, records, taxa_profiles, *dist_assignmet_para
                )

                write_fasta(
                    final_comp, records, basename, out_folder, min_taxa, out_suffix
                )
                s_i = "Done with fasta {}".format(fasta)
                s_d = "# retained components: {}".format(len(final_comp))
                if DEBUG:
                    print("Components are compatible, " + s_d)
                return (s_i, s_d)

            #  Components are incompatible.
            #  Return all components as separated alignments.
            #  Add unclustered seq using a distance method.
            else:
                final_comp = dist_assignment(
                    comp_1, records, taxa_profiles, *dist_assignmet_para
                )
                write_fasta(
                    final_comp, records, basename, out_folder, min_taxa, out_suffix
                )
                s_i = "Done with fasta {}".format(fasta)
                s_d = "# retained components: {}".format(len(final_comp))
                if DEBUG:
                    print("Components are incompatible, " + s_d)
                return (s_i, s_d)


