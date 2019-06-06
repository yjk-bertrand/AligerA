# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
Controle module for the AligerA phasing tool.
"""
__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import os
from multiprocessing import Manager

import aligera_scripts.phase_scripts as phase
from aligera_scripts.utilities import create_result_folder, process_future_fasta


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def get_starting_files(input_folder, output_folder, input_suffix, output_suffix):

    output_prefixes = [
        x.split(output_suffix)[0]
        for x in os.listdir(output_folder)
        if output_suffix in x
    ]
    starting_files = [
        x
        for x in os.listdir(input_folder)
        if input_suffix in x and x.split(input_suffix)[0] not in output_prefixes
    ]
    return starting_files


def PhaseAlleles(cfg, logger, cpu_number, **kargs):
    logger.info("Phasing alleles starting")
    in_folder = cfg["input_folder"]
    in_suffix = cfg["input_suffix"]
    out_suffix = cfg["output_suffix"]
    out_folder = cfg["output_folder"]
    create_result_folder(out_folder, logger)
    logger.debug("Checking the settings")
    temporary_folder = kargs["temp"]
    phase.check_settings(cfg, logger)
    starting = get_starting_files(in_folder, out_folder, in_suffix, out_suffix)

    starting_files = [x for x in starting if in_suffix in x]  # [:1]
    s = "There are {} input fasta files".format(len(starting_files))
    print(s)
    logger.info(s)
    samtools_v = phase.find_samtools_version(cfg)
    print("samtools version is: .X", samtools_v)

    if starting_files:
        manager = Manager()
        fastas = manager.Queue()
        result_dict = manager.dict()
        process_future_fasta(
            phase.phasing_contig,
            starting_files,
            result_dict,
            fastas,
            cpu_number,
            logger,
            cfg,
            temp=temporary_folder,
            samtools_version=samtools_v,
            tqdm_desc="phasing",
        )

    logger.info("Phasing alleles finished")


def main():
    # ===================================================================================
    #     for debugging
    # ===================================================================================
    pass


if __name__ == "__main__":
    main()

