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
import shutil
import sys
import shelve
import itertools
import collections

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from aligera_scripts.primary_scripts import (
    blast_xml_parser,
    filtering_components,
    find_dat_pairs,
    find_taxa,
    generate_final_dict,
    populate_shelve,
    retrievePrimaryComponents,
)
from aligera_scripts.utilities import run_subprocess, toPickle


# =======================================================================================
#                FUNCTIONS
# =======================================================================================

def get_large_files(fasta, high_limit):
    """
    Check that the number of sequences exceeds the higher limit
    """
    nbr_sequences = len(list(SeqIO.index(fasta, "fasta")))
    if nbr_sequences >= high_limit:
        return False
    return True


def get_align_size(fasta, low_limit, high_limit):
    """
    Check that the number of sequences belongs to the correct range
    """
    nbr_sequences = len(list(SeqIO.index(fasta, "fasta")))
    if nbr_sequences >= high_limit or nbr_sequences <= low_limit:
        return False
    return True


def filter_with_taxa(fasta, taxa_list, taxa_occurences, boolean):
    """
    Filter homology groups that have more occurrences of some taxa than 
    the preset limit. If boolean 'OR' is selected, one of the taxa in the list 
    that is above the limit suffice to identify the alignment. If 'AND' is 
    selected, all taxa in list need to be above the limit in order to identify 
    the group"""
    records = list(SeqIO.index(fasta, "fasta"))
    taxa = [record.split("|")[0] for record in records]
    if not taxa_list:
        taxa_list = taxa
    c = collections.Counter(taxa)
    if boolean == "OR":
        for t in taxa_list:
            if c[t] >= taxa_occurences:
                return False
        return True
    elif boolean == "AND":
        switch = False
        for t in taxa_list:
            if c[t] < taxa_occurences:
                switch = True
        return switch
    else:
        s = 'boolean parameter not recognized, must be either OR or AND'
        sys.exit(s)


def run_makeblastdb_second(sendbuf):
    """
    Runs makeblastdb
    """
    for fasta in sendbuf:
        blast_db_name = fasta.split(".")[0]
        cmd = "makeblastdb -in {fasta} -dbtype nucl -out {blast_db}".format(
            fasta=fasta, blast_db=blast_db_name
        )
        run_subprocess(cmd)


def fasta_to_blastdb_second(fasta_list, cfg):
    """
    Find all associations between a fasta file anda blast db
    """
    permutations = itertools.permutations(fasta_list, 2)
    associations = [(x[0], x[1].split(".fas")[0]) for x in permutations]
    return associations


def run_blast_second(sendbuf, cfg):
    """
    Run blast, different blast command from the primary scrips as it allows \
    for different blast programs
    """
    blast_prg = cfg["blast_program"]
    if blast_prg not in ["blastn", "tblastx"]:
        exception = "[error] Wrong choice of blast program: \
select either blastn or tblastx"
        raise Exception(exception)

    evalue = cfg["Evalue"]
    threads = cfg["num_threads"]
    max_targets = cfg["max_target_sequences"]
    for fasta_db in sendbuf:
        out_file = fasta_db[0].split(".fas")[0] + "_vs_" + fasta_db[1] + ".xml"
        blast_second(
            blast_prg, fasta_db[0], out_file, fasta_db[1], evalue, threads, max_targets
        )


def blast_second(blast_prg, seqs, out_file, db, evalue, threads, max_targets):
    """
    Run blast
    """
    cmd = "{blast_program} -outfmt 5 -db {database_name} -query {input_file} \
-out {output} -evalue {evalue} -num_threads {num_threads} \
-max_target_seqs {max_target_seqs}".format(
        blast_program=blast_prg,
        input_file=seqs,
        output=out_file,
        database_name=db,
        evalue=evalue,
        num_threads=threads,
        max_target_seqs=max_targets,
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        template = "An exception of type {0} occurred when trying to run command {2}. \
    Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)


def parse_pickle_second(sendbuf, min_align_length):
    """
    Pickles the outcome of xml parsing. The difference with the primary 
    script is that it is applied on a list of xmls.
    """
    for xml in sendbuf:
        if os.stat(xml).st_size != 0:
            blast_query, blast_database = (xml.strip(".xml")).split("_vs_")
            parsed_xml = blast_xml_parser(xml, blast_query, blast_database, min_align_length)
            pickle_name = blast_query + "_vs_" + blast_database + ".dat"
            toPickle(parsed_xml, pickle_name)


def parse_fasta_list_second(fasta_list):
    """
    Return a SeqIO index of a list of fasta files in a dictionary
    """
    fasta_dict = {}
    for fasta in fasta_list:
        try:
            name = fasta.split(".fas")[0]
            idx = SeqIO.index(fasta, "fasta")
            fasta_dict[name] = idx
        except:
            exception = "[Error] Problem parsing or spliting\
fasta file {}".format(
                fasta
            )
            raise Exception(exception)
    return fasta_dict


def fetch_sequences_second(final_dict, parsed_fasta_dict):
    """
    Assemble the sequences grouped into components into alignments. 
    Unlike the primary script, it does not write the alignment to disk
    """
    component_sequences = []
    for component in sorted(list(final_dict.keys())):
        sequencesList = []
        for taxon, contigs in final_dict[component].items():
            for contig in contigs:
                try:
                    contig_name = taxon + "|" + contig
                    record = parsed_fasta_dict[taxon][contig_name]
                    newRecord = SeqRecord(
                        seq=record.seq, name=record.name, id=record.name, description=""
                    )
                    sequencesList.append(newRecord)
                except:
                    exception = "[error] Problem fetching sequences \
                    with taxon: {0} contig:{1}".format(
                        taxon, contig
                    )
                    print(exception)
                    continue
        component_sequences.append(sequencesList)
    return component_sequences


def mk_secondary(fasta, cfg, **kargs):
    """
    Main function for SECONDARY ALIGNMENTS:
    Large alignments are disassembled and re-blasted with higher stringency
    """
    temp = kargs["temp"]
    remove_temp = kargs["remove_temp"]
    in_folder = cfg["input_folder"]
    in_suffix = cfg["input_suffix"]
    min_size = cfg["min_taxa_in_alignment"]
    max_size = cfg["limit_large_file"]
    filter_taxa = cfg["filtering_mode_taxa"]
    filter_size = cfg["filtering_mode_size"]
    passed_folder = cfg["output_folder_passed"]
    failed_folder = cfg["output_folder_failed"]
    filter_taxa_list = cfg["filtered_taxa_list"]
    filter_taxa_occurence = cfg["filtered_taxon_occurences"]
    filter_taxa_bool = cfg["filtered_taxon_boolean"]
    min_align_length = cfg["min_align_len"]
    print("Starting the work on fasta: {}".format(fasta))

    os.chdir(temp)
    fasta_temp_folder = fasta.split(in_suffix)[0]

    #  Create a folder for re-analysis in the temporary directory.
    if fasta_temp_folder not in [x for x in os.listdir(os.getcwd())]:
        os.mkdir(fasta_temp_folder)
    os.chdir(in_folder)
    shutil.copy(fasta, os.path.join(temp, fasta_temp_folder))
    #  The target fasta is not deleted in order not to break the integrity
    #   of datasets obtained at each step
    #    try:
    #        os.remove(fasta)
    #    except:
    #        exception = "[error] Unable to remove fasta: {}".format(fasta)
    #        raise Exception(exception)

    os.chdir(os.path.join(temp, fasta_temp_folder))
    try:
        records = list(SeqIO.parse(fasta, "fasta"))
        if len(records) == 0:
            exception = "[error] Fasta: {} is empty".format(fasta)
            raise Exception(exception)
    except:
        exception = "[error] Unable to read fasta: {}".format(fasta)
        raise Exception(exception)

    #  Create a dict holding taxa as key and all contigs records as values.
    #  Empty positions (in case of aligned sequences) are removed.
    taxa_dict = {}
    for record in records:
        taxon = record.name.split("|")[0]
        if taxon in list(taxa_dict.keys()):
            taxa_dict[taxon].append(record)
        else:
            taxa_dict[taxon] = [record]
    for taxon in taxa_dict:
        records = taxa_dict[taxon]
        new_records = []
        for record in records:
            new_record = SeqRecord(
                seq=Seq(str(record.seq).replace("-", ""), IUPAC.unambiguous_dna),
                name=record.name,
                id=record.id,
                description="",
            )
            new_records.append(new_record)
        SeqIO.write(new_records, taxon + ".fas", "fasta")

    temp_fastas = [x + ".fas" for x in taxa_dict.keys()]

    #  Sequences are sorted by taxa. New taxa sequences are turned into
    #  blastdb and blasted against each other.
    run_makeblastdb_second(temp_fastas)
    fasta_to_db = fasta_to_blastdb_second(temp_fastas, cfg)
    run_blast_second(fasta_to_db, cfg)

    #  Resulting xmls are parsed and components are extracted from reciprocal
    #   blast analyses.
    xml_files = [x for x in os.listdir(os.getcwd()) if ".xml" in x]
    if not xml_files:
        exception = "[error] No xml file to parse. An error occured during \
    blasting"
        raise Exception(exception)

    parse_pickle_second(xml_files, min_align_length)
    dat_files = [x for x in os.listdir(os.getcwd()) if ".dat" in x]
    taxa_pairs = find_dat_pairs(dat_files)
    taxa = find_taxa(dat_files)

    common_hits = shelve.open("common_hits_dict")
    for pair in taxa_pairs:
        populate_shelve(common_hits, pair)

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
    #  Remove temporary files
    for dat in dats:
        os.remove(dat)

    #  Extract sequences alignments that correspond to the components.
    #  Filter the resulting alignment according to either taxa or size criteria.
    parsed_fasta_dict = parse_fasta_list_second(temp_fastas)
    final_dict = shelve.open("final_dico_complete")
    filter_aligns, discard_aligns = filtering_components(final_dict, min_size)
    final_dict.close()
    comp_seqs = fetch_sequences_second(filter_aligns, parsed_fasta_dict)
    if not filter_taxa_list: filter_taxa_list = taxa

    for item in enumerate(comp_seqs):
        name_fasta = "alignment_{0}L{1}{suffix}".format(
            fasta_temp_folder.split("_")[1], item[0], suffix=cfg["output_suffix"]
        )
        SeqIO.write(sorted(item[1], key=lambda x: x.name), name_fasta, "fasta")

        if not filter_taxa and filter_size:
            if get_align_size(name_fasta, min_size, max_size + 1):
                try:
                    shutil.move(name_fasta, passed_folder)
                except:
                    print(
                        "[Error] Cannot move file {0} to folder {1}".format(
                            name_fasta, passed_folder
                        )
                    )
            else:
                try:
                    shutil.move(name_fasta, failed_folder)
                except:
                    print(
                        "[Error] Cannot move file {0} to folder {1}".format(
                            name_fasta, failed_folder
                        )
                    )

        elif not filter_size and filter_taxa:
            if filter_with_taxa(
                name_fasta, filter_taxa_list, filter_taxa_occurence, filter_taxa_bool
            ):
                
                try:
                    shutil.move(name_fasta, passed_folder)
                except:
                    print(
                        "[Error] Cannot move file {0} to folder {1}".format(
                            name_fasta, passed_folder
                        )
                    )                   
            else:
                try:
                    shutil.move(name_fasta, failed_folder)
                except:
                    print(
                        "[Error] Cannot move file {0} to folder {1}".format(
                            name_fasta, failed_folder
                        )
                    )
                        
        elif filter_size and filter_taxa:
            if filter_with_taxa(
                name_fasta, filter_taxa_list, filter_taxa_occurence, filter_taxa_bool
            ) and get_align_size(name_fasta, min_size, max_size + 1):
                
                try:
                    shutil.move(name_fasta, passed_folder)
                except:
                    print(
                        "[Error] Cannot move file {0} to folder {1}".format(
                            name_fasta, passed_folder
                        )
                    )                   
            else:
                try:
                    shutil.move(name_fasta, failed_folder)
                except:
                    print(
                        "[Error] Cannot move file {0} to folder {1}".format(
                            name_fasta, failed_folder
                        )
                    )    
        else:
            try:
                shutil.move(name_fasta, passed_folder)
            except:
                print(
                    "[Error] Cannot move file {0} to folder {1}".format(
                        name_fasta, passed_folder
                    )
                )            
            

    os.chdir(os.pardir)
    if remove_temp:
        shutil.rmtree(fasta_temp_folder, ignore_errors=True)
    os.chdir(in_folder)

    s_i = "Done retrieving secondary alignment from fasta: {}".format(fasta)
    s_d = "Fasta: {0} was split into {1} alignments".format(fasta, len(comp_seqs))

    return (s_i, s_d)
