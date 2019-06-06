# -*- coding: utf-8 -*-
#
#    Copyright (C) 2019 Yann J.K. BERTRAND: yjk_bertrand@ybertrand.org
#       All rights reserved.
#

"""
    Script design to parse a .ini file, check for field validity 
    type and/or presence) and return a parameters dictionary
"""
__version__ = "1.3"

# =======================================================================================
#               IMPORTS
# =======================================================================================

import sys
import os
import configparser
import ast
import datetime 
from aligera_scripts.utilities import create_folder, exiter


# =======================================================================================
#                FUNCTIONS
# =======================================================================================

cfg = configparser.ConfigParser(
    allow_no_value=True, interpolation=configparser.ExtendedInterpolation()
)
cfg_dict = {
    "GLOBAL": {},
    "PRIMARY": {},
    "PHASING": {},
    "SECONDARY": {},
    "STEP1": {},
    "STEP2": {},
    "STEP3": {},
    "STEP4": {},
    "STEP5": {},
    "STEP6": {},
    "STEP7": {},
}


def licence_validation(boolean):
    """
    Test whether the software licence is valid. 
    Continue if valid sys.exit() otherwise.
    """
    if boolean:
        present = datetime.date.today()
        validity = datetime.date(2019, 9, 1)
        if present > validity:
            s = """Your licence for AligerA has expired. Contact the project's developer
                to renew it:  yjk_bertrand@ybertrand.org"""
            exiter(s)


def select_analysis(analysis_dict):
    sections = list(analysis_dict.keys())
    "Check that at least one analysis has been selected"
    run_analysis = False
    for section in sections:
        variables = analysis_dict[section]
        for variable in variables:
            if cfg[section].getboolean(variable) is True:
                run_analysis = True
                break

    if not run_analysis:
        s = """[ERROR] No analysis selected .ini file {2} section [{0}]
        please set on an analysis in:
            {1}
        EXITING
        """.format(
            section, variables, ini_file
        )
        sys.exit(s)


def path_populate_cfg_dict(added_dict):
    "Test and add paths to cfg_dict"
    for section, variables in added_dict.items():
        for variable in variables:
            R = IsPath(cfg, section, variable)
            # print(R, variable)
            if R:
                cfg_dict[section][variable] = R
            else:
                s = """[ERROR] problem setting path in .ini file {3}
                in section: {0}
                parameter: {1}
                path: {2} cannot be located
                EXITING""".format(
                    section, variable, cfg[section][variable], ini_file
                )
                sys.exit(s)


def path_populate_no_check_cfg_dict(added_dict):
    "add paths to cfg_dict without checking their validity"
    for section, variables in added_dict.items():
        for variable in variables:
            cfg_dict[section][variable] = cfg[section][variable]


def boolean_populate_cfg_dict(added_dict):
    "Test and add booleans to cfg_dict"
    for section, variables in added_dict.items():
        for variable in variables:
            R = IsBool(cfg, section, variable)
            if R is not None:
                cfg_dict[section][variable] = R
            else:
                s = """[ERROR] ValueError in in .ini file {3}
                in section: {0}
                parameter: {1}
                variable: {2} is not a boolean
                set to [True, False, true, false, on, off, 0, 1]
                EXITING""".format(
                    section, variable, cfg[section][variable], ini_file
                )
                sys.exit(s)


def integer_populate_cfg_dict(added_dict):
    "Test and add integers to cfg_dict"
    for section, variables in added_dict.items():
        for variable in variables:
            R = IsInt(cfg, section, variable)
            if R is not None:
                cfg_dict[section][variable] = R
            else:
                s = """[ERROR] ValueError in .ini file {3}
                in section: {0}
                parameter: {1}
                variable: {2} is not an integer
                EXITING""".format(
                    section, variable, cfg[section][variable], ini_file
                )
                sys.exit(s)


def float_populate_cfg_dict(added_dict):
    "Test and add floats to cfg_dict"
    for section, variables in added_dict.items():
        for variable in variables:
            R = IsFloat(cfg, section, variable)
            if R is not None:
                cfg_dict[section][variable] = R
            else:
                s = """[ERROR] ValueError in .ini file {3}
                in section: {0}
                parameter: {1}
                variable: {2} is not a float
                EXITING""".format(
                    section, variable, cfg[section][variable], ini_file
                )
                sys.exit(s)


def list_populate_cfg_dict(added_dict):
    "add lists to cfg_dict"
    for section, variables in added_dict.items():
        for variable in variables:
            R = IsList(cfg, section, variable)
            if R:
                cfg_dict[section][variable] = R
            else:
                s = """[ERROR] ValueError in .ini file {3}
                in section: {0}
                parameter: {1}
                variable: {2} is not a list
                EXITING""".format(
                    section, variable, cfg[section][variable], ini_file
                )
                sys.exit(s)


def string_populate_cfg_dict(added_dict):
    "add strings to cfg_dict"
    for section, variables in added_dict.items():
        for variable in variables:
            cfg_dict[section][variable] = cfg[section][variable]


def IsFloat(cfg, section, variable):
    try:
        R = cfg[section].getfloat(variable)
        return R
    except ValueError:
        try:
            R = eval(cfg[section][variable])
            return R
        except:
            return None


def IsInt(cfg, section, variable):
    try:
        R = cfg[section].getint(variable)
        return R
    except:
        return None


def IsList(cfg, section, variable):
    try:
        R = ast.literal_eval(cfg[section][variable])
        if isinstance(R, list):
            return R
        else:
            return None
    except:
        return None


def IsBool(cfg, section, variable):
    try:
        R = cfg[section].getboolean(variable)
        return R
    except:
        return None


def IsPath(cfg, section, variable):
    if os.path.isdir(cfg[section][variable]):
        return cfg[section][variable]
    else:
        return None


def special_populate_cfg_dict(added_dict):
    "add predefined types"
    type_dict = {
        "path": IsPath,
        "float": IsFloat,
        "int": IsInt,
        "bool": IsBool,
        "list": IsList,
    }
    for section, variables in added_dict.items():
        for variable in variables:
            target = variable[0]
            list_types = variable[1]
            switch = False
            if cfg[section][target] == "None":
                cfg_dict[section][target] = None
                switch = True
            else:
                for t in list_types:
                    if t not in [type_dict.keys()]:
                        s = """[ERROR] ValueError in .ini file {3}
                        in section: {0}
                        parameter: {1}
                        variable: {2} is of unknown type
                        EXITING""".format(
                            section, variable[0], t, ini_file
                        )
                    func = type_dict[t]
                    R = func(cfg, section, target)
                    if R:
                        cfg_dict[section][target] = R
                        switch = True
                        break
            if switch is False:
                s = """[ERROR] ValueError in .ini file {3}
                in section: {0}
                parameter: {1} did not match any of
                the predefined types: {2}
                EXITING""".format(
                    section, target, list_types, ini_file
                )
                sys.exit(s)


def global_populate_cfg_dict(
    boolean_dict,
    analysis_dict,
    path_dict,
    path_dict_no_check,
    integer_dict,
    float_dict,
    list_dict,
    str_dict,
    special_dict={},
):
    "Global function for populating all types in cfg_dict"
    if boolean_dict:
        boolean_populate_cfg_dict(boolean_dict)
    if analysis_dict:
        select_analysis(analysis_dict)
    if path_dict:
        path_populate_cfg_dict(path_dict)
    if path_dict_no_check:
        path_populate_no_check_cfg_dict(path_dict_no_check)
    if integer_dict:
        integer_populate_cfg_dict(integer_dict)
    if float_dict:
        float_populate_cfg_dict(float_dict)
    if list_dict:
        list_populate_cfg_dict(list_dict)
    if str_dict:
        string_populate_cfg_dict(str_dict)
    if special_dict:
        special_populate_cfg_dict(special_dict)


def main(ini, script):
    licence_validation(True)
    cfg.read(ini)
    global ini_file
    ini_file = ini

    if not cfg.read(ini):
        s = """[ERROR] no .ini file named '{0}' in {1}
        EXITING""".format(
            ini_file, os.getcwd()
        )
        sys.exit(s)
    #  Add global parameters used throughout the programme
    global_dict = {"GLOBAL": ["aligera_folder"]}
    create_folder(cfg["GLOBAL"]["temporary_folder"])
    cfg_dict["GLOBAL"]["temporary_folder"] = cfg["GLOBAL"]["temporary_folder"]
    boolean_dict = {"GLOBAL": ["remove_temporary_files"]}
    path_populate_cfg_dict(global_dict)
    boolean_populate_cfg_dict(boolean_dict)
    # print([(x, y) for (x, y) in cfg_dict.items()])

    # Add local paremeters used within the tools
    if script == "PRIMARY":
        boolean_dict = {
            script: [
                "run_makeblastdb",
                "run_BLAST",
                "build_primary_groups",
                "keep_discarded_groups",
                "fasta_sanity_check",
                "parallel",
            ]
        }
        analysis_dict = {
            script: ["run_makeblastdb", "run_BLAST", "build_primary_groups"]
        }
        path_dict = {}
        path_dict_no_check = {
            script: {
                "input_folder": cfg[script]["input_folder"],
                "output_folder": cfg[script]["output_folder"],
                "blastn_databases_folder": cfg[script]["blastn_databases_folder"],
                "blastn_results_folder": cfg[script]["blastn_results_folder"],
                "discarded_groups_folder": cfg[script][
                    "discarded_groups_folder"
                ],
            }
        }
        integer_dict = {
            script: [
                "min_taxa_in_alignment",
                "num_threads",
                "max_target_sequences",
                "min_align_len",
            ]
        }
        float_dict = {script: ["Evalue"]}
        list_dict = {}
        str_dict = {script: ["input_suffix", "output_suffix"]}
        special_dict = {}
        global_populate_cfg_dict(
            boolean_dict,
            analysis_dict,
            path_dict,
            path_dict_no_check,
            integer_dict,
            float_dict,
            list_dict,
            str_dict,
            special_dict,
        )

    elif script == "PHASING":
        boolean_dict = {
            script: ["use_only_properly_paired_reads", "check_reads_validity"]
        }
        analysis_dict = {}
        path_dict = {script: ["input_folder", "bam_folder", "transcriptomes_folder"]}
        path_dict_no_check = {
            script: {
                "output_folder": cfg[script]["output_folder"],
                "samtools_path": cfg[script]["samtools_path"],
                "bcftools_path": cfg[script]["bcftools_path"],
                "vcfutils_path": cfg[script]["vcfutils_path"],
                "varscan_path": cfg[script]["varscan_path"],
            }
        }
        integer_dict = {
            script: [
                "maximum_sequences",
                "varscan_min_avg_qual",
                "min_coverage_for_snp",
            ]
        }
        float_dict = {script: ["varscan_min_fq", "validity_stingency"]}
        list_dict = {script: ["transcriptomes_with_paired_reads"]}
        str_dict = {script: ["input_suffix", "output_suffix"]}
        special_dict = {
            script: [
                ("taxa_excluded_from_phasing", ["list", "None"]),
                ("transcriptomes_with_paired_reads", ["list", "None"]),
            ]
        }
        global_populate_cfg_dict(
            boolean_dict,
            analysis_dict,
            path_dict,
            path_dict_no_check,
            integer_dict,
            float_dict,
            list_dict,
            str_dict,
            special_dict,
        )

    elif script == "SECONDARY":
        boolean_dict = {
            script: [
                "proceed_with_secondary_search",
                "filtering_mode_taxa",
                "filtering_mode_size",
                "fasta_sanity_check",
            ]
        }
        analysis_dict = {}
        path_dict = {script: ["input_folder"]}
        path_dict_no_check = {
            script: {
                "output_folder_passed": cfg[script]["output_folder_passed"],
                "output_folder_failed": cfg[script]["output_folder_failed"],
            }
        }
        integer_dict = {
            script: [
                "min_taxa_in_alignment",
                "limit_large_file",
                "filtered_taxon_occurences",
                "num_threads",
                "max_target_sequences",
                "min_align_len",
            ]
        }
        float_dict = {script: ["Evalue"]}
        list_dict = {}
        #list_dict = {script: ["filtered_taxa_list"]}
        str_dict = {
            script: [
                "input_suffix",
                "filtered_taxon_boolean",
                "blast_program",
                "output_suffix",
            ]
        }
        special_dict = {script: [("filtered_taxa_list", ["list"])]}
        global_populate_cfg_dict(
            boolean_dict,
            analysis_dict,
            path_dict,
            path_dict_no_check,
            integer_dict,
            float_dict,
            list_dict,
            str_dict,
            special_dict,
        )

    elif script == "ALIGERA PIPELINE":
        analysis_dict = {
            "STEP1": ["run_STEP1"],
            "STEP2": ["run_STEP2"],
            "STEP3": ["run_STEP3"],
            "STEP4": ["run_STEP4"],
            "STEP5": ["run_STEP5"],
            "STEP6": ["run_STEP6"],
            "STEP7": ["run_STEP7"],
        }

        select_analysis(analysis_dict)
        boolean_dict = {
            "STEP1": ["run_STEP1"],
            "STEP2": ["run_STEP2"],
            "STEP3": ["run_STEP3"],
            "STEP4": ["run_STEP4"],
            "STEP5": ["run_STEP5"],
            "STEP6": ["run_STEP6"],
            "STEP7": ["run_STEP7"],
        }

        boolean_populate_cfg_dict(boolean_dict)

        if cfg_dict["STEP1"]["run_STEP1"]:
            script = "STEP1"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {
                script: ["upper_sequence_limit", "MAFFT_upper_limit_addfragments"]
            }
            float_dict = {}
            list_dict = {}
            str_dict = {
                script: [
                    "input_suffix",
                    "output_suffix",
                    "MAFFT_parameters_small",
                    "MAFFT_parameters_large",
                ]
            }
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
            )

        if cfg_dict["STEP2"]["run_STEP2"]:
            script = "STEP2"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {script: ["window_length"]}
            float_dict = {
                script: ["ungap_proportion", "min_aa_proportion", "max_gap_proportion"]
            }
            list_dict = {}
            str_dict = {script: ["input_suffix", "output_suffix"]}
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
            )

        if cfg_dict["STEP3"]["run_STEP3"]:
            script = "STEP3"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {
                script: [
                    "kmer",
                    "max_Hamming_distance",
                    "consecutive_ambiguity",
                    "min_length_ambiguous_sequence",
                    "min_taxa_in_alignment",
                    "min_alignment_length",
                ]
            }
            float_dict = {script: ["max_ambiguity_proportion", "difference_proportion"]}
            list_dict = {script: ["outgroups"]}
            str_dict = {script: ["input_suffix", "output_suffix"]}
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
            )

        if cfg_dict["STEP4"]["run_STEP4"]:
            script = "STEP4"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {
                script: [
                    "min_sequences_overlap",
                    "min_taxa_in_alignment",
                    "max_alleles",
                    "dna_model",
                ]
            }
            float_dict = {
                script: [
                    "max_jaccard_value",
                    "min_association_ratio",
                    "max_distance",
                    "large_component_ratio",
                ]
            }
            list_dict = {script: ["kmers"]}
            str_dict = {script: ["input_suffix", "output_suffix"]}
            special_dict = {script: [("outgroups", ["list", "None"])]}
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
                special_dict,
            )

        if cfg_dict["STEP5"]["run_STEP5"]:
            script = "STEP5"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {
                script: [
                    "seq_overlap_splitting",
                    "min_taxa_in_alignment",
                    "seq_overlap_merging",
                    "dna_model",
                ]
            }
            float_dict = {
                script: [
                    "min_proportion_taxa_that_overlap",
                    "max_jaccard_value",
                    "paralog_allele_distance_ratio",
                    "large_component_ratio",
                    "max_alleles_distance",
                ]
            }
            list_dict = {}
            str_dict = {script: ["input_suffix", "input_format", "output_suffix"]}
            special_dict = {
                script: [
                    ("excluded_paralogs", ["list", "None"]),
                    ("outgroups", ["list", "None"]),
                ]
            }
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
                special_dict,
            )

        if cfg_dict["STEP6"]["run_STEP6"]:
            script = "STEP6"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {
                script: [
                    "lenght_small_fragment",
                    "length_flanking_gaps",
                    "window_length",
                    "length_tiny_fragments",
                ]
            }
            float_dict = {
                script: ["ungap_proportion", "min_aa_proportion", "max_gap_proportion"]
            }
            list_dict = {}
            str_dict = {script: ["input_suffix", "input_format", "output_suffix"]}
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
            )

        if cfg_dict["STEP7"]["run_STEP7"]:
            script = "STEP7"
            boolean_dict = {}
            analysis_dict = {}
            path_dict = {script: ["input_folder"]}
            path_dict_no_check = {script: ["output_folder"]}
            integer_dict = {
                script: ["kmer", "overlap_threshold", "mistmatch_cost", "dna_model"]
            }

            float_dict = {
                script: ["overlap_proportion", "max_alleles_distance", "overlap_cost"]
            }
            # list_dict = {script:["taxa_without_alleles"]}
            list_dict = {}
            str_dict = {script: ["input_suffix", "output_suffix", "input_format"]}
            special_dict = {script: [("taxa_without_alleles", ["list", "None"])]}
            global_populate_cfg_dict(
                boolean_dict,
                analysis_dict,
                path_dict,
                path_dict_no_check,
                integer_dict,
                float_dict,
                list_dict,
                str_dict,
                special_dict,
            )
    return cfg_dict


if __name__ == "__main__":
    script = "ALIGERA PIPELINE"
    cfg = main(script)
    print(cfg["STEP7"])
