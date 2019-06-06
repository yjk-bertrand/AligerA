# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
"""
Controling module for the  AligerA primary tool.
"""
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
import sys
import glob
from multiprocessing import Manager
import shelve

import aligera_scripts.primary_scripts as primaryScript
from aligera_scripts.utilities import (
    check_for_file,
    create_result_folder,
    is_fasta,
    process_future_fasta,
)

# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def runMakeblastdb(cfg, logger, cpu_number):
    logger.info("Running blastdb")
    in_folder = cfg["input_folder"]
    in_suffix = cfg["input_suffix"]
    blastn_db_folder = cfg["blastn_databases_folder"]
    create_result_folder(blastn_db_folder, logger)
    os.chdir(in_folder)
    finished_files = [x.split(".nin")[0] for x in os.listdir(blastn_db_folder)]
    starting_files = sorted(
        [
            x
            for x in os.listdir(os.getcwd())
            if in_suffix in x and x.split(in_suffix)[0] not in finished_files
        ]
    )
    if cfg["fasta_sanity_check"]:
        for fasta in starting_files:
            is_fasta(fasta)
    logger.info("There are {} input fasta files".format(len(starting_files)))
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            primaryScript.run_makeblastdb,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Running blastdb",
        )

    logger.info("run_makeblastdb finished")


def runBLAST(cfg, logger, cpu_number):
    logger.info("Running BLAST")
    in_folder = cfg["input_folder"]
    blastn_dbs = cfg["blastn_databases_folder"]
    in_suffix = cfg["input_suffix"]
    out_folder = cfg["blastn_results_folder"]
    create_result_folder(out_folder, logger)
    if not os.path.exists(blastn_dbs) or not os.listdir(blastn_dbs):
        exception = "[error]: either blastdb folder {} is empty or \
it has not been created by running 'run_makeblastdb'".format(
            blastn_dbs
        )
        raise Exception(exception)
        sys.exit()
        
    if not os.path.exists(in_folder) or not os.listdir(in_folder):
        exception = "[error]: either fasta input folder {} is empty or \
it does not exist".format(
            in_folder
        )
        raise Exception(exception)
        sys.exit()
    os.chdir(in_folder)
    starting_files = sorted([x for x in os.listdir(os.getcwd()) if in_suffix in x])

    if cfg["fasta_sanity_check"]:
        for fasta in starting_files:
            primaryScript.is_fasta(fasta)

    logger.info("There are {} input fasta files".format(len(starting_files)))
    fasta_to_db = primaryScript.fasta_to_blastdb(starting_files, cfg)

    logger.info("Starting BLAST...")
    if fasta_to_db:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            primaryScript.run_blastn,
            fasta_to_db,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            tqdm_desc="Running BLAST",
        )

    logger.info("run_BLAST finished")


def parseXMLs(cfg, logger, cpu_number, temporary_folder):
    logger.info("Parsing xmls")

    in_folder = cfg["blastn_results_folder"]
    if not os.path.exists(in_folder) or not [
        x for x in os.listdir(in_folder) if ".xml" in x
    ]:
        exception = "[error]: either blastn_results_folder {} is empty or \
it has not been created by running 'run_BLAST'".format(
            in_folder
        )
        raise Exception(exception)
        sys.exit()
    os.chdir(in_folder)
    starting_files = sorted([x for x in os.listdir(os.getcwd()) if ".xml" in x])

    logger.info("There are {} input xml files".format(len(starting_files)))
    logger.debug("Parsing xml files")
    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            primaryScript.parse_pickle,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            temporary_folder,
            min_align_length=cfg["min_align_len"],
            tqdm_desc="Parsing xmls",
        )
    os.chdir(temporary_folder)
    dat_files = [x for x in os.listdir(os.getcwd()) if ".dat" in x]
    taxa_pairs = primaryScript.find_dat_pairs(dat_files)
    common_hits = shelve.open("common_hits_dict")
    for pair in taxa_pairs:
        primaryScript.populate_shelve(common_hits, pair)
    common_hits.close()
    logger.info("Done parsing xmls")


def componentsSequential(cfg, logger, temporary_folder):
    logger.info("Building components sequentially")
    os.chdir(temporary_folder)
    check_for_file("common_hits_dict.dat")
    primaryScript.getComponentsSequential(cfg, logger)
    # ===================================================================================
    for shl in glob.glob("common_hits_dict.*"):
        os.remove(shl)
    # ===================================================================================
    logger.info("Done with computing components")


def componentsParallel(cfg, logger, cpu_number, temporary_folder):
    logger.info("Building components in parallel")
    os.chdir(temporary_folder)
    check_for_file("common_hits_dict.dat")
    primaryScript.getComponentsParallel(cfg, logger, cpu_number)
    # ===================================================================================
    for shl in glob.glob("common_hits_dict.*"):
        os.remove(shl)
    # ===================================================================================
    logger.info("Done with computing components")


def makePrimaryAssortments(cfg, logger, temporary_folder):
    logger.info("Retrieving the primary sequence assortments")
    out_folder = cfg["output_folder"]
    create_result_folder(out_folder, logger)
    discard_aligns = cfg["discarded_groups_folder"]
    create_result_folder(discard_aligns, logger)
    primaryScript.mk_primary(cfg, logger, temporary_folder)


if __name__ == "__main__":
    pass
