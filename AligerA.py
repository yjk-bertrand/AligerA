# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

__version__ = "1.3"


# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import os
import time
import datetime
import argparse
import logging
import logging.handlers
from multiprocessing import cpu_count

import aligera_scripts.config_parser as config_parser
import aligera_scripts.utilities as utilities

# =======================================================================================
#               COMMAND LINE PARSING
# =======================================================================================


def create_parser():
    parser = argparse.ArgumentParser(formatter_class=argparse.MetavarTypeHelpFormatter)
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    parser.add_argument(
        "--tool",
        "-t",
        type=str,
        choices=["primary", "secondary", "phasing", "pipeline"],
        help="Select a tool among the AligerA's applications.",
        required=True,
    )
    parser.add_argument(
        "--verbose",
        "-v",
        help="Increase logging verbosity from INFO to DEBUG.",
        action="store_true",
    )
    parser.add_argument(
        "--cpu",
        type=int,
        help="Set the number of processors for \
                        parallelization.\n\
                        Default: all available processors",
        default=cpu_count(),
    )
    parser.add_argument(
        "--log",
        "-l",
        type=str,
        help="""Set the logging file.\n\
                        Default: 'AligerA_log.txt'""",
        default="AligerA_log.txt",
    )
    parser.add_argument(
        "--ini",
        "-i",
        type=str,
        help="""Set the .ini file.\n\
                        Default: 'aligera.ini'""",
        default="aligera.ini",
    )
    return parser


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def ping():
    args_tool = "pipeline"
    # args_tool = "phasing"
    args_verbose = True
    args_cpu = 1
    args_log = "aligera_pipeline_log.txt"
    args_ini = "aligera.ini"
    return args_tool, args_verbose, args_cpu, args_log, args_ini


# =======================================================================================
#                MAIN
# =======================================================================================


def main():
    # ===================================================================================
    #     set logger and parse arguments
    # ===================================================================================
    start_time = time.time()
    parser = create_parser()
    args = parser.parse_args()

    if args.verbose:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO
    LOG_FILENAME = args.log
    logger = logging.getLogger("MyLogger")

    handler = logging.handlers.RotatingFileHandler(
        LOG_FILENAME, maxBytes=2000000, backupCount=5
    )
    logger.addHandler(handler)
    start = "Pipeline started at date & time {0}".format(time.strftime("%c"))
    print(start)
    logger.setLevel(logging_level)
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging_level)

    available_cpus = cpu_count()
    if args.cpu:
        cpu_number = args.cpu
        if cpu_number > available_cpus:
            sys.exit(
                "[Error] Exiting:\
            more CPU requested than the available number of {} CPUS".format(
                    available_cpus
                )
            )
        if cpu_number == 0:
            sys.exit("[Error] Exiting: '--cpu 0' is not a useful parameter")
    print("cpu_number = ", cpu_number)
    print("verbosity", logging_level)

    # ===================================================================================
    #      get global parameters
    # ===================================================================================

    cfg = config_parser.main(args.ini, "GLOBAL")
    global aligera_folder
    aligera_folder = cfg["GLOBAL"]["aligera_folder"]
    global temporary_folder
    temporary_folder = cfg["GLOBAL"]["temporary_folder"]
    utilities.create_result_folder(temporary_folder, logger)
    global clean_temp_files
    clean_temp_files = cfg["GLOBAL"]["remove_temporary_files"]
    if clean_temp_files:
        logger.info("Deleting files in temporary folder")
        utilities.cleanTemporaryFolder(temporary_folder)
        os.chdir(aligera_folder)
    if args.verbose: 
        for (k, v) in cfg["GLOBAL"].items():
            print(k,":",v)
    logger.debug(cfg["GLOBAL"])
    
    # ===================================================================================
    #         run primary
    # ===================================================================================

    if args.tool == "primary":
        import aligera_scripts.AligerA_primary as primary  

        s = "\nSTARTING PRIMARY TOOL"
        print(s)
        logger.info(s)
        cfg = config_parser.main(args.ini, "PRIMARY")
        if cfg["PRIMARY"]["run_makeblastdb"]:
            s = "\nBuilding BLAST db"
            print(s)
            logger.info(s)
            primary.runMakeblastdb(cfg["PRIMARY"], logger, cpu_number)

        if cfg["PRIMARY"]["run_BLAST"]:
            s = "\nStarting BLAST"
            print(s)
            logger.info(s)
            primary.runBLAST(cfg["PRIMARY"], logger, cpu_number)

        if cfg["PRIMARY"]["build_primary_groups"]:
            s = "\nParsing BLAST xmls"
            print(s)
            logger.info(s)
            primary.parseXMLs(
                cfg["PRIMARY"], logger, cpu_number, temporary_folder
            )
            if not cfg["PRIMARY"]["parallel"]:
                s = "\nBuilding primary groups sequentially"
                print(s)
                logger.info(s)
                primary.componentsSequential(
                cfg["PRIMARY"], logger, temporary_folder
                )
            else:
                s = "\nBuilding primary groups in parallel"
                print(s)
                logger.info(s)
                primary.componentsParallel(
                cfg["PRIMARY"], logger, cpu_number, temporary_folder
                )
                
            s = "\nCreating groups"
            print(s)
            logger.info(s)
            primary.makePrimaryAssortments(
                cfg["PRIMARY"], logger, temporary_folder
            )
                
        if cfg["PRIMARY"]["build_primary_groups"] and clean_temp_files:
            logger.info("Deleting files in temporary folder")
            utilities.cleanTemporaryFolder(temporary_folder)

    # ===================================================================================
    #         run secondary
    # ===================================================================================

    if args.tool == "secondary":
        import aligera_scripts.AligerA_secondary as secondary  

        s = "\nSTARTING SECONDARY TOOL"
        print(s)
        logger.info(s)
        cfg = config_parser.main(args.ini, "SECONDARY")
        secondary.MakeSecondaryAssortments(
            cfg["SECONDARY"], logger, cpu_number, temp=temporary_folder, 
            remove_temp=clean_temp_files,
        )
        if (
            cfg["SECONDARY"]["proceed_with_secondary_search"]
            and clean_temp_files
        ):
            logger.info("Deleting files in temporary folder")
            utilities.cleanTemporaryFolder(temporary_folder)

    # ===================================================================================
    #         run phasing
    # ===================================================================================

    if args.tool == "phasing":
        import aligera_scripts.AligerA_phase as phase  

        s = "\nSTARTING PHASING"
        print(s)
        logger.info(s)
        cfg = config_parser.main(args.ini, "PHASING")
        phase.PhaseAlleles(cfg["PHASING"], logger, cpu_number, temp=temporary_folder)

        if clean_temp_files:
            logger.info("Deleting files in temporary folder")
            utilities.cleanTemporaryFolder(temporary_folder)

    # ===================================================================================
    #         run pipeline
    # ===================================================================================

    if args.tool == "pipeline":
        import aligera_scripts.AligerA_pipeline as pipeline  

        s = "\nSTARTING PIPELINE"
        print(s)
        logger.info(s)
        cfg = config_parser.main(args.ini, "ALIGERA PIPELINE")

        if cfg["STEP1"]["run_STEP1"]:
            s = "\nSTARTING STEP1"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP1(cfg["STEP1"], logger, cpu_number)

        if cfg["STEP2"]["run_STEP2"]:
            s = "\nSTARTING STEP2"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP2(cfg["STEP2"], logger, cpu_number)

        if cfg["STEP3"]["run_STEP3"]:
            s = "\nSTARTING STEP3"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP3(cfg["STEP3"], logger, cpu_number)

        if cfg["STEP4"]["run_STEP4"]:
            s = "\nSTARTING STEP4"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP4(cfg["STEP4"], logger, cpu_number)

        if cfg["STEP5"]["run_STEP5"]:
            s = "\nSTARTING STEP5"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP5(cfg["STEP5"], logger, cpu_number)

        if cfg["STEP6"]["run_STEP6"]:
            s = "\nSTARTING STEP6"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP6(cfg["STEP6"], logger, cpu_number)

        if cfg["STEP7"]["run_STEP7"]:
            s = "\nSTARTING STEP7"
            print(s)
            logger.info(s)
            pipeline.pipeline_STEP7(cfg["STEP7"], logger, cpu_number)

    secs = time.time() - start_time
    logger.info("Pipeline completed  in {} hours:minutes:seconds".format(
            datetime.timedelta(seconds=int(secs))))

    logger.removeHandler(handler)
    logging.shutdown()


if __name__ == "__main__":
    print("\nSTARTING ANALYSES")
    main()
    print("ANALYSES COMPLETED\n")
