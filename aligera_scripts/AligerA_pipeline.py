# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Controle module for the AligerA pipeline tool.
"""
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
from multiprocessing import Manager

from Bio import SeqIO

import aligera_scripts.STEP1_scripts as STEP1
import aligera_scripts.STEP2_scripts as STEP2
import aligera_scripts.STEP3_scripts as STEP3
import aligera_scripts.STEP4_scripts as STEP4
import aligera_scripts.STEP5_scripts as STEP5
import aligera_scripts.STEP6_scripts as STEP6
import aligera_scripts.STEP7_scripts as STEP7
from aligera_scripts.utilities import (
    create_result_folder,
    is_fasta,
    process_future_fasta
)


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def find_n_seqs(fasta, low_limit, high_limit):
    """
    Calculate the number of sequences in a fasta file
    """
    nbr_sequences = len(list(SeqIO.index(fasta, "fasta")))
    if nbr_sequences >= high_limit or nbr_sequences <= low_limit:
        return False
    else:
        return True


def find_n_taxa(fasta, low_limit):
    """
    Calculate the number of taxa that don't have empty sequences in a fasta file.
    """
    records = SeqIO.index(fasta, "fasta")
    seq_names = list(records.keys())
    taxa = set([name.split("|")[0] for name in seq_names if len(records[name].seq) > 0])
    if len(taxa) >= low_limit:
        return True
    else:
        return False


def cleaning_dnadist(logger):
    """
    Remove all temporary files associated with the distance calculations.
    """
    to_clean = [x for x in os.listdir(os.getcwd()) if "dnadist" in x]
    logger.info("Cleaning 'dnadist' files")
    for f in to_clean:
        os.remove(f)


def pipeline_STEP1(cfg, logger, cpu_number):
    s = "Working with STEP 1: running MAFFT"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    output_folder = cfg["output_folder"]
    out_suffix = cfg["output_suffix"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    unwanted_files = [x for x in os.listdir(os.getcwd()) if "_temp_aligned.fasta" in x]

    # Remove files from a previously aborted run .
    for f in unwanted_files:
        os.remove(f)
    finished_files = [x.split(out_suffix)[0] for x in os.listdir(output_folder)]
    starting_files = sorted(
        [
            x
            for x in os.listdir(os.getcwd())
            if cfg["input_suffix"] in x
            and find_n_seqs(x, 1, cfg["upper_sequence_limit"])
            and "_core" not in x
            and "_addit" not in x
            and x.split(out_suffix)[0] not in finished_files
        ]
    )

    for fasta in starting_files:
        is_fasta(fasta)

    # Size threshold for switching from small_fastas alignment to large_fastas
    size_threshold = cfg["MAFFT_upper_limit_addfragments"]
    small_fastas = [x for x in starting_files if find_n_seqs(x, 1, size_threshold + 1)]
    logger.debug("there are {} small_fastas".format(len(small_fastas)))

    large_fastas = [
        x for x in starting_files if not find_n_seqs(x, 1, size_threshold + 1)
    ]
    logger.debug("there are {} large_fastas".format(len(large_fastas)))

    if small_fastas:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        logger.info("Running Mafft on small fasta")
        process_future_fasta(
            STEP1.run_MAFFT_small,
            small_fastas,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Mafft on small files",
        )
    if large_fastas:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        logger.info("Running Mafft on large fasta")
        process_future_fasta(
            STEP1.run_MAFFT_large,
            large_fastas,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Mafft on large files",
        )
    logger.info("STEP 1 finished")


def pipeline_STEP2(cfg, logger, cpu_number):
    s = "Working with STEP 2: ALIGNMENT TRIMMING"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    output_folder = cfg["output_folder"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    starting_files = get_starting_files(
        input_folder, output_folder, in_suffix, out_suffix
    )
    logger.debug("There are {} input fasta files".format(len(starting_files)))
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            STEP2.trim,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Trimming",
        )
    logger.info("STEP 2 finished")


def pipeline_STEP3(cfg, logger, cpu_number):
    s = "Working with STEP 3: ALIGNMENT CLEANING"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    output_folder = cfg["output_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    starting_files = get_starting_files(
        input_folder, output_folder, in_suffix, out_suffix
    )
    logger.debug("There are {} input fasta files".format(len(starting_files)))
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            STEP3.find_unspliced_segments,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Cleaning",
        )
    logger.info("STEP 3 finished")


def pipeline_STEP4(cfg, logger, cpu_number):
    s = "Working with STEP 4: OUTPARALOGS SEPARATION"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    output_folder = cfg["output_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    min_taxa = cfg["min_taxa_in_alignment"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    starting = get_starting_files(input_folder, output_folder, in_suffix, out_suffix)
    starting_files = [
        x for x in starting if in_suffix in x and find_n_taxa(x, min_taxa)
    ]
    logger.info("There are {} input fasta files".format(len(starting_files)))

 
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            STEP4.outparalog_separation,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Outparalogs",
        )
      
    cleaning_dnadist(logger)
    logger.info("STEP 4 finished")


def pipeline_STEP5(cfg, logger, cpu_number):
    s = "Working with STEP 5: INPARALOGS SEPARATION"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    output_folder = cfg["output_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    in_format = cfg["input_format"]
    min_taxa = cfg["min_taxa_in_alignment"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    starting = get_starting_files(input_folder, output_folder, in_format, out_suffix)
    starting_files = [
        x for x in starting if in_suffix in x and find_n_taxa(x, min_taxa)
    ]
    logger.info("There are {} input fasta files".format(len(starting_files)))
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            STEP5.inparalog_separation,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Inparalogs",
        )
    cleaning_dnadist(logger)
    logger.info("STEP 5 finished")


def pipeline_STEP6(cfg, logger, cpu_number):
    s = "Working with STEP 6: SECOND ALIGNMENT CLEANING"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    output_folder = cfg["output_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    in_format = cfg["input_format"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    starting = get_starting_files(input_folder, output_folder, in_format, out_suffix)
    starting_files = [x for x in starting if in_suffix in x]
    logger.info("There are {} input fasta files".format(len(starting_files)))
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            STEP6.trim_clean,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Cleaning",
        )
    logger.info("STEP 6 finished")


def pipeline_STEP7(cfg, logger, cpu_number):
    s = "Working with STEP 7: Building final sequences"
    print(s)
    logger.info(s)
    input_folder = cfg["input_folder"]
    output_folder = cfg["output_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    in_format = cfg["input_format"]
    create_result_folder(output_folder, logger)
    os.chdir(input_folder)
    starting = get_starting_files(input_folder, output_folder, in_format, out_suffix)
    starting_files = [x for x in starting if in_suffix in x]
    logger.info("There are {} input fasta files".format(len(starting_files)))
    """
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            STEP7.assemble_alleles,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="STEP 7",
        )
        
    """
    for fasta in starting_files:
        STEP7.assemble_alleles(fasta, cfg)
 
    cleaning_dnadist(logger)
    logger.info("STEP 7 finished")


def get_starting_files(input_folder, output_folder, input_suffix, output_suffix):
    output_prefixes = [
        x.split(output_suffix)[0]
        for x in os.listdir(output_folder)
        if output_suffix in x
    ]
    starting_files = [
        x
        for x in os.listdir(os.getcwd())
        if input_suffix in x and x.split(input_suffix)[0] not in output_prefixes
    ]
    return starting_files


def main():
    # ===================================================================================
    #     for testing
    # ===================================================================================
    pass


if __name__ == "__main__":
    main()
