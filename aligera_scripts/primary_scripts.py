#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
"""
Changes with previous version:
    Components extraction has been modified: previous versions relied on
    creating an undirected graph and computing connected components using
    DFS in networkx. Starting with v.1.1 component computation is performed
    using the Union-Find algorithm with path compression and Union-by-rank [1].
    Running time is O(n.alpha(m, n)) where alpha is the very slowy growing
    inverse Ackerman function [2].
    Is introduced a parallel version with p processors that should run in
    O(n.alpha(m, n)/p) + O(mlog2p) |V| = m, |E| = n.
   
    [1] Cormen, T.H., Leiserson, C.E., Rivest, R.L., Stein, C. 2009: 
    Introduction to Algorithms (3 ed.). Cambridge, Massachusetts: MIT Press. 
    [2] Fredrik Manne , Md. Mostofa Ali Patwary, A scalable parallel union-find 
    algorithm for distributed memory computers, Proceedings of the 8th international 
    conference on Parallel processing and applied mathematics: Part I,
    September 13-16, 2009, Wroclaw, Poland 

"""
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
import sys
import pickle
import shelve
import itertools
from multiprocessing import Manager


from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from aligera_scripts.utilities import (
    grouper,
    is_fasta,
    process_future_shl,
    toPickle,
    run_subprocess,
)

from aligera_scripts.union_find import UnionFind
from aligera_scripts.merge_components import merge

# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def fasta_to_blastdb(fasta_list, cfg):
    """
    Find all associations between a fasta file and a blast db
    """
    blast_db_folder = cfg["blastn_databases_folder"]
    in_suffix = cfg["input_suffix"]
    blast_dbs = [x.split(".nin")[0] for x in os.listdir(blast_db_folder)]
    fasta_prefix = [x.split(in_suffix)[0] for x in fasta_list]
    missing_dbs = [x + in_suffix for x in fasta_prefix if x not in blast_dbs]
    if missing_dbs:
        s = "[error] the following fasta file do not have a corresponding \
blastdb: {}".format(
            ", ".join(missing_dbs)
        )
        sys.exit(s)
    perm = itertools.permutations(fasta_list, 2)
    associations = [
        (x[0], os.path.join(blast_db_folder, x[1].split(in_suffix)[0])) for x in perm
    ]
    return associations


def run_makeblastdb(fasta, cfg, **kargs):
    """
    Runs makeblastdb
    """
    blastn_databases_folder = cfg["blastn_databases_folder"]
    blast_db_name = os.path.join(blastn_databases_folder, fasta.split(".")[0])
    cmd = "makeblastdb -in {fasta} -dbtype nucl -out {blast_db_name}".format(
        fasta=fasta, blast_db_name=blast_db_name
    )
    run_subprocess(cmd)
    s_i = "Done building blastdb for fasta {}".format(fasta)
    s_d = None
    return (s_i, s_d)


def run_blastn(fasta_db, cfg, **kargs):
    """
    Runs blastn on pair of  fasta: fasta_db[0] against blastdb: fasta_db[1]
    """
    blastn_results = cfg["blastn_results_folder"]
    in_suffix = cfg["input_suffix"]
    evalue = cfg["Evalue"]
    num_threads = cfg["num_threads"]
    max_target_seqs = cfg["max_target_sequences"]
    fasta = fasta_db[0]
    blastdb = fasta_db[1]
    out_file = os.path.join(
        blastn_results,
        fasta.split(in_suffix)[0] + "_vs_" + os.path.basename(blastdb) + ".xml",
    )
    blastn(fasta, out_file, blastdb, evalue, num_threads, max_target_seqs)
    s_i = None
    s_d = "Done blasting fasta: {0} against blastdb: {1}".format(fasta, blastdb)

    return (s_i, s_d)


def blastn(seqFile, out_file, db, evalue, num_threads, max_target_seqs):
    blastResultFile = out_file
    cmd = "blastn -outfmt 5 -db {database_name} -query {input_file} \
-out {output} -evalue {evalue} -num_threads {num_threads} \
-max_target_seqs {max_target_seqs}".format(
        input_file=seqFile,
        output=blastResultFile,
        database_name=db,
        evalue=evalue,
        num_threads=num_threads,
        max_target_seqs=max_target_seqs,
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        template = "An exception of type {0} occurred when trying to run command {2}. \
    Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)


def blast_xml_parser(xml, blast_query, blast_db, min_align_length):
    """
    Parse the xml from running blast. Extract all hits on target
    """
    L_multiple_hits = []
    for block in list(NCBIXML.parse(open(xml))):
        if block.alignments:
            query = block.query
            hitsList = []
            for alignment in block.alignments:
                if [h for h in alignment.hsps if h.align_length > min_align_length]:
                    hitsList.append(alignment.hit_def)
            L_multiple_hits.extend(
                [
                    tuple(sorted([blast_query + "|" + query, blast_db + "|" + x]))
                    for x in hitsList
                ]
            )
    return L_multiple_hits


def parse_pickle(xml, temporary_folder, min_align_length=300, **kargs):
    """
    Pickles the outcome of xml parsing
    """
    if os.stat(xml).st_size != 0:
        blast_query, blast_database = (xml.split(".xml")[0]).split("_vs_")
        parsed_xml = blast_xml_parser(
            xml, blast_query, blast_database, min_align_length
        )
        pickle_name = blast_query + "_vs_" + blast_database + ".dat"
        toPickle(parsed_xml, os.path.join(temporary_folder, pickle_name))
    s_i = None
    s_d = "Done parsing xml: {0}".format(xml)
    return (s_i, s_d)


def find_dat_pairs(dat_files):
    """
    find reciprocal blast i.e. A_vs_B.xml is associated with B_vs_A.xml
    """
    pairs = []
    stack = dat_files[1:]
    target = dat_files[0]
    while stack:
        for item in stack:
            if (
                set(target[:-4].split("_vs_")).intersection(
                    set(item[:-4].split("_vs_"))
                )
                == set(target[:-4].split("_vs_"))
                and len(set(target[:-4].split("_vs_"))) == 2
            ):
                pairs.append((target, item))
        target = stack[0]
        stack = stack[1:]
    return pairs


def find_taxa(dat_files):
    """
    return the identity of taxa present in the blast output
    """
    taxa = []
    for x in dat_files:
        taxa.extend(x[:-4].split("_vs_"))
    return list(set(taxa))


def populate_shelve(shelve, pair):
    """
    Find the reciprocal hits
    """
    i, j = pair
    L1 = pickle.load(open(i, "rb"))
    L1_sorted = [tuple(sorted(x)) for x in L1]
    L2 = pickle.load(open(j, "rb"))
    L2_sorted = [tuple(sorted(x)) for x in L2]
    name_split = sorted(i[:-4].split("_vs_"))
    name = "{0}|{1}".format(name_split[0], name_split[1])
    common = set(L1_sorted).intersection(set(L2_sorted))
    shelve[name] = common


def retrievePrimaryComponents(k, dict_elts, shl, **kargs):
    """
    Compute components with Union-Find 
    from pairs in open shelve (shl) values,
    k partitions of shl keys,
    dict_elts is a dictionary with partitioin number as key and
    lists of shelve keys (species) as values
    """
    keys = dict_elts[k]
    with UnionFind() as UV:
        for i in keys:
            for item in shl[i]:
                UV.union(item)
    return k, UV.getComponents()


def getComponentsSequential(cfg, logger):
    """
    Compute connected components sequentially using Union-Find. 
    Save the components in a shelve
    """
    print("Computing components sequentially")
    try:
        common_hits = shelve.open("common_hits_dict")
    except:
        s = "[Error] cannot open shelve 'common_hits_dict' in {}".format(os.getcwd())
        sys.exit(s)
    dict_elts = {0: list(common_hits.keys())}
    components = retrievePrimaryComponents(0, dict_elts, common_hits)
    common_hits.close()
    primary_comp_dict = dict([(i, j) for i, j in enumerate(components[1])])
    generate_final_dict(primary_comp_dict)
    dats = [
        x
        for x in os.listdir(os.getcwd())
        if ".dat" in x
        and "_vs_" in x
        and "common_hits_dict" not in x
        and "final_dico_complete" not in x
    ]

    for dat in dats:
        os.remove(dat)


def getComponentsParallel(cfg, logger, cpu_number):
    """
    Compute spanning forests for p partitions of input data with Union-Find
    in parallel on p processors.
    Merge the forest pair-wise sequentially to obtain all connected components.
    Save the components in a shelve
    """
    print("Computing components in parallel")
    try:
        common_hits = shelve.open("common_hits_dict")
    except:
        s = "[Error] cannot open shelve 'common_hits_dict' in {}".format(os.getcwd())
        sys.exit(s)
    cpu_number = min(cpu_number, len(list(common_hits.keys())))
    if cpu_number == 1:
        logger.debug("Only one pair of taxa: switching to getComponentsSequential()")
        common_hits.close()
        getComponentsSequential(cfg, logger)
    else:
        sendbuf = grouper(cpu_number, list(common_hits.keys()))
        manager = Manager()
        queue = manager.Queue()
        result_dict = manager.dict()
        process_future_shl(
            retrievePrimaryComponents,
            sendbuf,
            result_dict,
            queue,
            cpu_number,
            common_hits,
            tqdm_desc="Components",
        )
        common_hits.close()
        forests = list(result_dict.values())
        components = merge(forests)
        primary_comp_dict = dict([(i, j) for i, j in enumerate(components)])
        generate_final_dict(primary_comp_dict)
        dats = [
            x
            for x in os.listdir(os.getcwd())
            if ".dat" in x
            and "_vs_" in x
            and "common_hits_dict" not in x
            and "final_dico_complete" not in x
        ]
        for dat in dats:
            os.remove(dat)


def generate_final_dict(component_dict):
    """
    Generate the dictionary that assign a name (key)
    to each component (value)
    """
    final_dict = shelve.open("final_dico_complete")
    i = 0
    for component in component_dict.values():
        keys = list(set([x.split("|")[0] for x in component]))
        empty_lists = [[] for x in range(len(keys))]
        name = "alignment_{0}".format(i)
        D = dict(zip(keys, empty_lists))
        for seq in component:
            D[seq.split("|")[0]] = D[seq.split("|")[0]] + [
                seq.split("|")[-1].split()[0]
            ]
        final_dict[name] = D
        i += 1
    final_dict.close()


def filtering_components(final_dict, min_taxa_in_alignment):
    """
    Filter component based on number of taxa
    """
    print("Filtering components...")
    filtered_alignments_keys = []
    for name, sequences in final_dict.items():
        if len(list(sequences.keys())) > min_taxa_in_alignment:
            filtered_alignments_keys.append(name)

    filtered_alignments = dict(
        zip(
            ["alignment_{}".format(i) for i in range(len(filtered_alignments_keys))],
            [final_dict[x] for x in filtered_alignments_keys],
        )
    )
    discarded_keys = [x for x in final_dict.keys() if x not in filtered_alignments_keys]
    discarded_alignments = dict(
        zip(
            ["alignment_{}".format(i) for i in range(len(discarded_keys))],
            [final_dict[x] for x in discarded_keys],
        )
    )
    print("Done filtering components...")
    return filtered_alignments, discarded_alignments


def parse_fasta_list(fastas, cfg):
    """
    Return a SeqIO index of a list of fasta files in a dictionary
    """
    fasta_dict = {}
    in_suffix = cfg["input_suffix"]
    for fasta in fastas:
        try:
            name = fasta.split(in_suffix)[0]
            idx = SeqIO.index(fasta, "fasta")
            fasta_dict[name] = idx
        except:
            exception = "[Error] Problem parsing or spliting \
fasta file {}".format(
                fasta
            )
            raise Exception(exception)
    return fasta_dict


def fetch_sequences(final_dict, parsed_fasta_dict, cfg):
    """
    Assemble the sequences grouped into components into alignments
    """
    out_suffix = cfg["output_suffix"]
    for component in sorted(list(final_dict.keys())):
        sequencesList = []
        for taxon, contigs in final_dict[component].items():
            for contig in contigs:
                try:
                    record = parsed_fasta_dict[taxon][contig]
                    new_seq = str(record.seq).replace("N", "").replace("n", "")
                    newRecord = SeqRecord(
                        seq=Seq(new_seq, alphabet=IUPAC.ambiguous_dna),
                        name=taxon + "|" + record.name,
                        id=taxon + "|" + record.name,
                        description="",
                    )
                    sequencesList.append(newRecord)
                except:
                    exception = "problem with taxon: {0} contig: {1}".format(
                        taxon, contig
                    )
                    raise Exception(exception)
                    sys.exit()
        SeqIO.write(
            sorted(sequencesList, key=lambda x: x.name), component + out_suffix, "fasta"
        )


def mk_primary(cfg, logger, temporary_folder):
    """
    Fetch the sequences that correspond to the component and \
    build unaligned primary assorted sequences files
    """
    in_folder = cfg["input_folder"]
    in_suffix = cfg["input_suffix"]
    min_size = cfg["min_taxa_in_alignment"]

    if not os.path.exists(in_folder):
        s = "[error] input folder {} is either empty or does not exist".format(
            in_folder
        )
        sys.exit(s)
    fasta_files = [x for x in os.listdir(in_folder) if in_suffix in x]

    if not fasta_files:
        s = "[error] input folder {} is either empty or does not exist".format(
            in_folder
        )
        sys.exit(s)

    os.chdir(in_folder)

    if cfg["fasta_sanity_check"]:
        for fasta in fasta_files:
            is_fasta(fasta)

    parsed_fasta_dict = parse_fasta_list(fasta_files, cfg)
    logger.info("Done parsing fasta")
    final_dict = shelve.open(os.path.join(temporary_folder, "final_dico_complete"))
    filter_aligns, discard_aligns = filtering_components(final_dict, min_size)
    logger.info("Number of groups before filtering: {}".format(len(final_dict)))
    logger.info("Number of groups after filtering: {}".format(len(filter_aligns)))
    os.chdir(cfg["output_folder"])
    logger.info("Fetching sequences...")
    fetch_sequences(filter_aligns, parsed_fasta_dict, cfg)

    if cfg["keep_discarded_groups"]:
        os.chdir(cfg["discarded_groups_folder"])
        fetch_sequences(discard_aligns, parsed_fasta_dict, cfg)

    logger.info("Done fetching sequences")
    final_dict.close()
    logger.info("Done generating primary groups")


def errors(msg):
    print(msg)
    sys.exit(msg)
