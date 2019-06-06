# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#
"""
Controle module for the  AligerA secondary tool.
"""
__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
import sys
import shutil
import re
from multiprocessing import Manager

from aligera_scripts.utilities import create_result_folder, process_future_fasta
import aligera_scripts.secondary_scripts as secondary
import aligera_scripts.primary_scripts as primary

# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def MakeSecondaryAssortments(cfg, logger, cpu_number, **kargs):
    logger.info("Retrieving the Secondary sequence assortments")
    in_folder = cfg["input_folder"]
    in_suffix = cfg["input_suffix"]
    max_size = cfg["limit_large_file"]
    filter_taxa = cfg["filtered_taxa_list"]
    filter_taxa_occurences = cfg["filtered_taxon_occurences"]
    filter_taxa_bool = cfg["filtered_taxon_boolean"]
    out_passed = cfg["output_folder_passed"]
    create_result_folder(out_passed, logger)
    out_failed = cfg["output_folder_failed"]
    create_result_folder(out_failed, logger)
    temporary_folder = kargs["temp"]
    remove_temp = kargs["remove_temp"]

    out_files = [re.split("[L][0-9]+", x)[0] for x in os.listdir(out_passed)] + [
        re.split("[L][0-9]+", x)[0] for x in os.listdir(out_failed)
    ]
    os.chdir(in_folder)
    starting_files = sorted(
        [
            x
            for x in os.listdir(os.getcwd())
            if (
                in_suffix in x
                and "_core" not in x
                and "_addit" not in x
                and "_temp_" not in x
                and os.stat(x).st_size != 0
                and x not in out_files
            )
        ]
    )

    if cfg["fasta_sanity_check"]:
        logger.debug("Performing sanity check")
        for fasta in starting_files:
            primary.is_fasta(fasta)

    logger.debug("Searching for large fasta files")
    starting_fasta = "fasta files:" + ",".join(starting_files)
    logger.debug(starting_fasta)
    filtering_modes = "filtering_mode_size: {0}, filtering_mode_taxa: {1}".format(
            cfg["filtering_mode_size"], cfg["filtering_mode_taxa"])
    logger.debug(filtering_modes)
    if cfg["filtering_mode_size"] and not cfg["filtering_mode_taxa"]:
        large_fastas = [
            x
            for x in starting_files
            if not secondary.get_large_files(x, max_size + 1)
        ]
    elif cfg["filtering_mode_taxa"] and not cfg["filtering_mode_size"]:
        large_fastas = [
            x
            for x in starting_files
            if not secondary.filter_with_taxa(
                x, filter_taxa, filter_taxa_occurences, filter_taxa_bool
            )
        ]
    elif cfg["filtering_mode_taxa"] and cfg["filtering_mode_size"]:
        large_fastas = [
            x
            for x in starting_files
            if (not secondary.filter_with_taxa(
                x, filter_taxa, filter_taxa_occurences, filter_taxa_bool)
                and not secondary.get_large_files(x, max_size + 1)
            )
        ]         
    else:
        large_fastas = []
        
    small_fastas = list(set(starting_files).difference(set(large_fastas)))
    for x in small_fastas:
        try:
            shutil.copy(x, out_passed)
        except:
            print("[Error]: unable to move file {0} to folder {1}".format(x, out_passed))
            
    logger.info("There are {} small fasta files that have been validated".format(len(small_fastas)))
    logger.info("There are {} large fasta files that need to be processed".format(len(large_fastas)))
    

    if cfg["proceed_with_secondary_search"] and large_fastas:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            secondary.mk_secondary,
            large_fastas,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            temp=temporary_folder,
            tqdm_desc="2nd search",
            remove_temp=remove_temp,
        )

    logger.debug("Done MakeSecondaryAssortments")


def main():
    # ===================================================================================
    #     for testing
    # ===================================================================================
    pass


if __name__ == "__main__":
    main()
