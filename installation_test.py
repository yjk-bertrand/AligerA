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
import sys
import subprocess
import logging
import logging.handlers


def check_version(cmd_str, tool_name, version, logger):
    """
    Throw exception if run command fails
    return stdout and stderr
    """
    v = str(version)
    process = subprocess.Popen(
        cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    stdout, stderr = process.communicate()
    
    logger.info("\nTesting {} installation".format(tool_name))

    if not stdout:
        logger.info("{} is NOT properly installed".format(tool_name))
        logger.info(stderr)
    else:
        tool_version = stdout.decode().split()[1]
        if not tool_version.startswith(v):
            logger.info("Wrong {tool} version: must be {v}.X".format(tool=tool_name,
                        v=v))
        else:
            logger.info("{} installation.......ok".format(tool_name))
            return False
    return True


def run_subprocess(cmd_str):
    """
    return stdout and stderr
    """
    process = subprocess.Popen(
        cmd_str, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    stdout_str, stderr_str = process.communicate()
    return stdout_str, stderr_str


def main(ini_file):
    print("Running installation test with ini file: {}".format(ini_file))
    LOG_FILENAME = "installation_report.txt"
    logger = logging.getLogger("MyLogger")
    logging_level = logging.INFO
    handler = logging.handlers.RotatingFileHandler(
        LOG_FILENAME, maxBytes=2000000, backupCount=5
    )
    logger.addHandler(handler)
    logger.setLevel(logging_level)
    stderrhandler = logging.StreamHandler(sys.stderr)
    stderrhandler.setFormatter(logging.Formatter("    %(message)s"))
    stderrhandler.setLevel(logging_level)
    logger.info("STARTING installation tests")
    try:
        import aligera_scripts.config_parser as config_parser
    except:
        logger.info("Could not load module config_parse.py \n \
AligerA scritps are missing or wrong python version...exiting")
        sys.exit("[Error] Cannot perform installation test, see log.")
    
    #  Tests
    switch = False
    logger.info("\nTesting python installation")
    python_out, python_err = run_subprocess("python -V")
    if not python_err.decode().startswith("Python 3"):
        logger.info("Python is NOT properly installed, requires Python 3.x")
        switch = True
    else:
        logger.info("Python installation.......ok")
    
    switch = check_version("conda -V", "anaconda", 4, logger)
    
    switch = check_version("blastn -version", "BLAST", 2, logger)
    
    logger.info("\nTesting samtools installation")
    cfg = config_parser.main(ini_file, "PHASING")
    samtools_str = cfg["PHASING"]["samtools_path"]
    samtools_out, samtools_err = run_subprocess(samtools_str)
    if not samtools_err.decode().startswith("\nProgram: samtools"):
        logger.info("Samtools is NOT properly installed")
        logger.info(samtools_err.decode())
        switch = True
    else:
        logger.info("Samtools installation.......ok")
    
    bcftools_str = cfg["PHASING"]["bcftools_path"] + " --version"
    switch = check_version(bcftools_str, "bcftools", 1, logger)
   
    logger.info("\nTesting vcfutils installation")
    vcfutils_str = cfg["PHASING"]["vcfutils_path"]
    vcfutils_out, vcfutils_err = run_subprocess(vcfutils_str)
    if not vcfutils_err.decode().startswith("\nUsage:   vcfutils.pl"):
        logger.info("vcfutils is NOT properly installed")
        logger.info(vcfutils_err.decode())
        switch = True
    else:
        logger.info("vcfutils installation.......ok")
        
    logger.info("\nTesting mafft installation")       
    mafft_out, mafft_err = run_subprocess("mafft --version") 
    if not mafft_err.decode().startswith("v7."):
        logger.info("mafft is NOT properly installed, needs version 7.X")
        logger.info(mafft_err.decode())
        switch = True    
    else:
        logger.info("mafft installation.......ok")
        
    logger.info("\nTesting java installation")        
    java_out, java_err = run_subprocess("java -version")   
    java_switch = False
    if not java_err.decode().startswith("java version"):
        logger.info("java is NOT properly installed, needs version 7.X")
        logger.info(java_err.decode())
        logger.info("CANNOT test for Varscan installation")
        switch = True  
        java_switch = True
    else:
        logger.info("java installation.......ok")
        
    if not java_switch:
        varscan_str = cfg["PHASING"]["varscan_path"]
        logger.info("\nTesting VarScan installation")
        varscan_out, varscan_err = run_subprocess("java -jar {}".format(varscan_str))   
        if not varscan_err.decode().startswith("VarScan v2"):
            logger.info("VarScan is NOT properly installed, needs version 2.X")
            logger.info(varscan_err.decode())
            switch = True  
        else:
            logger.info("VarScan installation.......ok")   
 
    logger.info("\nTesting emboss installation") 
    emboss_out, emboss_err = run_subprocess("distmat -help")
    if not emboss_err.decode().split()[10].startswith("EMBOSS:6"):
        logger.info("EMBOSS is NOT properly installed, needs version 6.X")
        logger.info(emboss_err.decode())
        switch = True
    else:
        logger.info("EMBOSS installation.......ok")  
        
    if switch:
        print("ERROR in third party installation, check log file for details")
    else:
        print("Third party installation.......ok")
    logger.removeHandler(handler)
    logging.shutdown()

if __name__ == "__main__":
    if "installation_report.txt" in os.listdir("."):
        os.remove("installation_report.txt")
    if not len(sys.argv) == 2:
        sys.exit("please provide an INI file as argument")
    ini_file = sys.argv[1]
    main(ini_file)

