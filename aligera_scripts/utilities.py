#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Utility functions
    isfloat(value)

    create_result_folder(folder, logger)
    
    create_folder(folder)
    
    cleanTemporaryFolder(temporary_folder)
    
    grouper(n, iterable)
    
    run_subprocess(cmd_str)
    
    process_future_shl(script_to_run, list_elts, result_dict, queue, cpu_number, shl)
                                                   
    process_future_fasta(script_to_run, fasta_list, result_dict,queue, cpu_number, 
                         logger, cfg) 
    group_tasks(queue, result_dict)
    
    check_for_file(file)
    
    toPickle(data, file_name)
    
    is_fasta(fasta)
    
    get_alleles(taxa, outgroups)
    
    sequence_profiler(record)
    
    get_overlap_score(taxa, taxa_profiles)
    
    remove_ambiguities(seq, char_to_remove)
    
    calculate_distance_matrix(fasta, cfg)
    
    distance_matrix_parser(distance_matrix)
    
"""
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
import sys
import shutil
import subprocess
import pickle
from multiprocessing import Queue, Manager
import concurrent.futures

from tqdm import tqdm

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def exiter(message):
    sys.exit(message)

def isfloat(value):
    """
    Check for float type
    """
    try:
        float(value)
        return True
    except ValueError:
        return False


def create_result_folder(folder, logger):
    """
    Create folders if not present
    """
    if not os.path.exists(folder):
        if logger:
            logger.info("Could not find folder {}...try to create it".format(folder))
        try:
            os.mkdir(folder)
            if logger:
                logger.info("Folder creation successful")
        except:
            if logger:
                logger.info("[ERROR] Failed to create folder {}\n".format(folder))


def create_folder(folder):
    """
    Create folders if not present
    """
    if not os.path.exists(folder):
        try:
            os.mkdir(folder)
        except:
            s = "[ERROR] Failed to create folder {}\n".format(folder)
            sys.exit(s)


def cleanTemporaryFolder(temporary_folder):
    os.chdir(temporary_folder)
    for x in os.listdir(os.getcwd()):
        try:
            os.remove(x)
        except:
            shutil.rmtree(x)


def grouper(n, iterable):
    return [iterable[i::n] for i in range(n)]


def run_subprocess(cmd_str, get_stdout=False, get_stderr=False):
    """
    Throw exception if run command fails
    return either stdout or stderr
    """
    process = subprocess.Popen(
        cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    stdout_str, stderr_str = process.communicate()
    if process.returncode != 0:
        raise Exception(
            "Failed to run '%s'\n%s%sNon-zero exit status %s"
            % (cmd_str, stdout_str, stderr_str, process.returncode)
        )
        sys.exit()
    if get_stdout:
        return stdout_str
    if get_stderr:
        return stderr_str
    return


def process_future_shl(
    script_to_run, list_elts, result_dict, queue, cpu_number, shl, **kargs
):
    """
    compute asynchronously with multiprocessing and ProcessPoolExecutor
    list_elts is a list of partitions of shelve keys
    result_dict is the dictionary linked to the Manager
    shl is an open shelve.
    """
    description = "Anonymous Task"
    if kargs and "tqdm_desc" in list(kargs.keys()):
        description = kargs["tqdm_desc"]
    dict_elts = dict(enumerate(list_elts))
    keys = list(dict_elts.keys())
    for k in keys:
        queue.put(k)

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=cpu_number
    ) as group_processes:
        for i in range(queue.qsize()):
            group_processes.submit(group_tasks, queue, result_dict)

    with concurrent.futures.ProcessPoolExecutor(max_workers=cpu_number) as processes:
        future_tasks = {
            processes.submit(script_to_run, k, dict_elts, shl): k
            for k in result_dict.keys()
        }
        kwargs = {
            "total": len(future_tasks),
            "unit": "nap",
            "unit_scale": False,
            "leave": True,
            "desc": description,
        }

        for future in tqdm(concurrent.futures.as_completed(future_tasks), **kwargs):
            """
            if future.exception() is not None:
                s = "[Error] fasta file {0} has generated an exception: \
                {1}".format(file_name, future.exception())                
                sys.exit(s)
            """
            result_dict[future.result()[0]] = future.result()[1]


def process_future_fasta(
    script_to_run, fasta_list, result_dict, queue, cpu_number, logger, cfg, **kargs
):
    """
    process the fasta files using multiprocessing and ProcessPoolExecutor
    """
    description = "Anonymous Task"
    if kargs and "tqdm_desc" in list(kargs.keys()):
        description = kargs["tqdm_desc"]

    for fasta in fasta_list:
        queue.put(fasta)
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=cpu_number
    ) as group_fasta_processes:
        for i in range(queue.qsize()):
            group_fasta_processes.submit(group_tasks, queue, result_dict)

    with concurrent.futures.ProcessPoolExecutor(
        max_workers=cpu_number
    ) as fasta_processes:
        future_tasks = {
            fasta_processes.submit(script_to_run, fasta, cfg, **kargs): fasta
            for fasta in result_dict.keys()
        }
        kwargs = {
            "total": len(future_tasks),
            "unit": "nap",
            "unit_scale": False,
            "leave": True,
            "desc": description,
        }

        for future in tqdm(concurrent.futures.as_completed(future_tasks), **kwargs):
            file_name = future_tasks[future][0]
            
            if future.exception() is not None:
                s = "[Error] fasta file named: {0} has generated an exception: \
                {1}".format(
                    file_name, future.exception()
                )
                sys.exit(s)
            else:
                if future.result()[0]:
                    logger.info(future.result()[0])
                if future.result()[1]:
                    logger.debug(future.result()[1])


def group_tasks(queue, result_dict):
    """
    Populate queue,
    Add keys to result_dict
    """
    try:
        q = queue.get(True, 0.05)
        result_dict[q] = None
    except queue.Empty:
        print("Nothing to be done, queue is empty")


def check_for_file(file):
    """
    Check that file exists and is not empty.
    """
    if file not in os.listdir(os.getcwd()) or os.stat(file).st_size == 0:
        exception = "[error]: The file {0} is either \
empty or it could not be found in the folder: {1}".format(
            file, os.getcwd()
        )
        raise Exception(exception)
        sys.exit()


def toPickle(data, file_name):
    """
    pickle: dump data in a file
    """
    with open(file_name, "wb") as f:
        pickle.dump(data, f)


def is_fasta(fasta):
    """
    Check that the file is a non empty fasta file
    """
    try:
        L = list(SeqIO.parse(fasta, "fasta"))
        if len(L) == 0:
            s = "[error] Problem with fasta {}: the file is empty".format(fasta)
            raise Exception(s)
            sys.exit()
    except:
        s = "[error] Fail to open file {}, make sure the file is in fasta \
format".format(
            fasta
        )
        raise Exception(s)
        sys.exit()
    return True


def get_alleles(taxa, outgroups):
    """
    Identify the alleles.
    The alleles are identified if they correspond to
    the allele_0/allele_1 designation obtained from the phasing step.
    taxa and outgroups are lists of sequence names.
    """
    taxa_with_alleles = list(
        set(
            [
                x[0]
                for x in [
                    y.split("allele") for y in taxa if y.split("|")[0] not in outgroups
                ]
                if x[0] + "allele_0" in taxa and x[0] + "allele_1" in taxa
            ]
        )
    )

    return [(p + "allele_0", p + "allele_1") for p in taxa_with_alleles]


def sequence_profiler(record):
    """
    Turn each sequence into bit array i.e. a string of 0 and 1, 0 for "-", 
    1 else. Aimed to provide a fast and easy way to calculate sequence overlap.
    Returns a dict with record name as key and bitarray as value.
    interp. record (dict)   A dictionary with taxa names as keys and as values their
                            aligned sequences in str format.
    """
    taxa_profiles = {
        k: list(map(lambda x: 0 if x == "-" else 1, list(v))) for k, v in record.items()
    }
    return taxa_profiles


def get_overlap_score(taxa, taxa_profiles):
    """
    Compute sequence overlap from sequence profile.
    interp.   taxa (list)    a pair of keys of the taxa_profiles dictionary
              taxa_profiles (dict)    A dictionary that holds sequence names as keys
                                      and as values their sequences transformed into
                                      bitarrays with '0' for missing/ambiguous base
                                      '1' otherwise.
    """
    item_overlap = sum(
        a == b for a, b in zip(taxa_profiles[taxa[0]], taxa_profiles[taxa[1]]) if a == 1
    )
    return item_overlap


def remove_ambiguities(seq, char_to_remove):
    """
    Remove unwanted characters from string.
    'char_to_remove' is a str of unwanted char
    """
    mapper = dict.fromkeys(i for i in range(sys.maxunicode) if chr(i) in char_to_remove)
    return seq.translate(mapper)


def calculate_distance_matrix(fasta, cfg):
    """
    Use EMBOSS distmat program to compute the distance matrix
    between sequences in an alignment. A temporary fasta file without ambiguity (Nn-) is
    created (fasta+'_for_dnadist.fasta').
    interp.
        fasta (Bio.SeqIO)   Sequences object.
        cfg (dict)          Dictionary with parameters
    """
    try:
        records = list(SeqIO.parse(fasta, "fasta"))
    except:
        exception = "[Error]: Cannot open alignment {}\n\
        Verify that the file is in fasta format".format(
            fasta
        )
        return exception
    basename = fasta.split(cfg["input_suffix"])[0]
    new_records = []
    for record in records:
        cleaned_seq = remove_ambiguities(str(record.seq), "Nn-")
        if len(cleaned_seq) > 100:
            new_records.append(
                SeqRecord(
                    seq=record.seq,
                    name="_WITH_".join((record.name).split("|")),
                    id="_WITH_".join((record.name).split("|")),
                    description="",
                )
            )
    SeqIO.write(new_records, basename + "_for_dnadist.fasta", "fasta")
    cmd = "distmat -sequence {outfile} -nucmethod {model}\
           -outfile {outfile}.distmat".format(
        model=cfg["dna_model"], outfile=basename + "_for_dnadist.fasta"
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        template = "An exception of type {0} occurred when trying to run command {2}. \
    Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)
    os.remove(basename + "_for_dnadist.fasta")


def distance_matrix_parser(distance_matrix):
    """
    Parser that returns a list of pair-wise distances
    from the distance matrix
    """
    try:
        with open(distance_matrix, "r") as f:
            lines = list(f.readlines()[8:])
    except:
        s = "[Error] Unable to open distance matrix {} \n\
        Exiting".format(
            distance_matrix
        )
        raise Exception(s)
        sys.exit()
    splitted_lines = [x.split("\t") for x in lines]
    name_index = {}
    for index, item in enumerate(splitted_lines):
        name = item[-1].split(" ")[0].replace("_WITH_", "|")
        name_index[index] = name
    distances = []
    for taxon in splitted_lines[:]:
        target_taxon = taxon[-1].split(" ")[0].replace("_WITH_", "|")
        for index, distance in enumerate(taxon[1:-2]):
            if isfloat(distance) and distance not in ["-nan", "", "nan"]:
                if float(distance) != 0.0:
                    distances.append(
                        (target_taxon, name_index[index], {"distance": float(distance)})
                    )
    try:
        os.remove(distance_matrix)
    except:
        s = "[Error] Unable to remove distance matrix file {} \n\
        Exiting".format(
            distance_matrix
        )
        raise Exception(s)
    return distances
