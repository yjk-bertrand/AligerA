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
import subprocess
import shelve

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

from aligera_scripts.utilities import run_subprocess

import pysam


# =======================================================================================
#                CLASSES
# =======================================================================================


class RawRead(object):
    """
    Pair or inpair read.
    """

    def __init__(self, value):
        parsed = value.split("\t")
        self.name = parsed[0]
        self.length = int(parsed[8])
        self.sequence = parsed[9]
        self.start = int(parsed[3])
        self.end = self.start + len(self.sequence)

    def __getitem__(self, index):
        if index not in range(self.start, self.end):
            raise IndexError("list index out of range")
        new_index = index - self.start
        return self.sequence[new_index]


class Read(RawRead):
    """
    Find reads that have been paired.
    """

    def __init__(self, value):
        if len(value) == 1:
            RawRead.__init__(self, value[0])
        else:
            self.read1 = RawRead(value[0])
            self.read2 = RawRead(value[1])
            self.name = self.read1.name
            self.start = min(self.read1.start, self.read2.start)
            self.end = max(self.read1.end, self.read2.end)
            self.length = self.end - self.start

            if self.start == self.read1.start:
                self.sequence = (
                    self.read1.sequence
                    + "-" * (self.read2.start - self.read1.end)
                    + self.read2.sequence
                )
            else:
                self.sequence = (
                    self.read2.sequence
                    + "-" * (self.read1.start - self.read2.end)
                    + self.read1.sequence
                )


# =======================================================================================
#                FUNCTIONS
# =======================================================================================


def find_sequence_number(fasta, cfg):
    """
    Check whether the sequence number in the alignment is below a preset threshold.
    """
    nbr_sequences = len(list(SeqIO.index(fasta, "fasta")))
    if nbr_sequences >= cfg["maximum_sequences"]:
        return False
    return True


def sanity_checker(cfg, logger):
    """
    Check that all bam files that correspond to the transcriptomes are present.
    """
    in_suffix = cfg["input_suffix"]
    fasta_origin = cfg["transcriptomes_folder"]
    if cfg["taxa_excluded_from_phasing"] is None:
        cfg["taxa_excluded_from_phasing"] = []
    #   Names of fastas derived from transcriptome sequencing
    trans_bams = [
        x.split(in_suffix)[0]
        for x in os.listdir(fasta_origin)
        if in_suffix in x and x not in cfg["taxa_excluded_from_phasing"]
    ]

    bams = [
        x.split(".bam")[0].split("_sorted")[0]
        for x in os.listdir(cfg["bam_folder"])
        if ".bam" in x and "bam.bai" not in x
    ]

    if trans_bams and bams:
        if set(trans_bams).intersection(set(bams)) == set(trans_bams) and set(
            trans_bams
        ).intersection(set(bams)) == set(bams):
            logger.info("all bam files are present")
        else:
            missing_bams = list(
                set(trans_bams).difference(set(trans_bams).intersection(set(bams)))
            )
            error = "Warning there are some missing bam files: {}".format(missing_bams)
            logger.error(error, exc_info=True)
            raise Exception(error)
    else:
        error = "[error] Either or both the bam or the transcriptome folders are empty"
        logger.error(error, exc_info=True)
        raise Exception(error)


def find_samtools_version(cfg):
    """
    Determine whether the samtools version is 0.X or 1.X.
    """
    process = subprocess.Popen(
        cfg["samtools_path"],
        stdout=subprocess.PIPE,
        shell=False,
        stderr=subprocess.PIPE,
    )
    out, err = process.communicate()
    version = None
    for line in err.decode().split("\n"):
        if line.startswith("Version"):
            version = int(line.split("Version: ")[1][0])
            break
    if version is None or version not in [0, 1]:
        sys.exit("Unable to recognize the SAMtools version...Exiting")

    return version


def calculate_mean_coverage_SD(bam_file, cfg):
    """
    On a bam file compute the mean coverage and the standard deviation 
    which provide for each contig the read coverage.
    """
    print("mean coverage not implemented")
    pass


def extract_reference(record, cfg):
    """
    Takes the fasta files with all contigs and extract the right sequence.
    """
    taxon = record.name.split("|")[0]
    contig = record.name.split("|")[1]
    new_record = SeqRecord(record.seq, name=contig, id=contig, description="")
    SeqIO.write(new_record, taxon + "_WITH_" + contig + "_ref.fas", "fasta")
    taxon_name = taxon + "_WITH_" + contig + "_ref.fas"
    cmd = "{samtools} faidx {}".format(taxon_name, samtools=cfg["samtools_path"])
    run_subprocess(cmd)


def call_snp(record, cfg):
    """
    Call snps using mpileup and varscan.
    """
    taxon = record.name.split("|")[0]
    contig = record.name.split("|")[1]

    cmd = "{samtools} mpileup -f {0}_WITH_{1}_ref.fas {0}_WITH_{1}_sorted.bam \
> {0}_WITH_{1}_ref.mpileup".format(
        taxon, contig, samtools=cfg["samtools_path"]
    )

    run_subprocess(cmd)

    #  Case 1. mpileup has returned an empty file. Do not run Varscan.
    if os.stat("{0}_WITH_{1}_ref.mpileup".format(taxon, contig)).st_size == 0:
        # print("mpileup empty")
        f = open("{0}_WITH_{1}_spns.txt".format(taxon, contig), "wb")
        f.close()
        return
    #  Case 2. the file from mpileup is not empty. Run Varscan.
    else:
        cmd = "java -jar {0} pileup2snp {1}_WITH_{2}_ref.mpileup \
--min-avg-qual 30 --min-var-freq {varscan_min_fq} --min-coverage \
{coverage} > {1}_WITH_{2}_spns.txt".format(
            cfg["varscan_path"],
            taxon,
            contig,
            coverage=cfg["min_coverage_for_snp"],
            varscan_min_fq=cfg["varscan_min_fq"],
        )
        try:
            run_subprocess(cmd)
        except:
            exception = "Unable to run varscan with taxon: {taxon}, \
contig: {contig} with command: {command}".format(
                taxon=taxon, contig=contig, command=cmd,
            )
            raise Exception(exception)
        return


def find_chromosome_name(sam_file, samtools):
    """
    Find the name of the reference sequence used for the mapping.
    """
    cmd = [samtools, "idxstats", sam_file]
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    sam_name = process.communicate()[0].decode().split("\t")[0]
    return sam_name


def database_sam(sam_file, sam_name, sam_start, sam_end, shelve_name):
    """
    Index the reads contained in the sam file by creating read objects. 
    Store the objects in a shelve used as a database.
    """
    samfile = pysam.AlignmentFile(sam_file, "rb")
    reads_list = []
    for read in samfile.fetch(sam_name, sam_start, sam_end):
        reads_list.append(str(read))
    reads_dict = {}
    for r in reads_list:
        read_name = r.split("\t")[0]
        if read_name not in list(reads_dict.keys()):
            reads_dict[read_name] = [r]
        else:
            reads_dict[read_name].append(r)
    reads_objects_list = []
    for item in reads_dict.values():
        reads_objects_list.append(Read(item))

    db = shelve.open(shelve_name)
    for read in reads_objects_list:
        db[read.name] = read
    db.close()


def count_overlapping_reads(snp_pairs, db):
    """
    Count the number of reads that span two adjacent snps.
    """
    snp_pairs_found = []
    remaining_snps = snp_pairs

    for read in db.values():
        for pair in remaining_snps:
            if read.start <= pair[0] and read.end > pair[1]:
                if read[pair[0]] != "-" and read[pair[1]] != "-":
                    snp_pairs_found.append(pair)
        remaining_snps = list(set(snp_pairs).difference(set(snp_pairs_found)))
    return len(snp_pairs_found)


def check_validity_stringency(samtools, bam_file, stringency, reference, snp_report):
    """
    Check whether the snps contained in pair-reads sufficiently overlap the total snps. 
    A stringency of 0 requires no overlap and 1 requires a 100% overlap. 
    snp_report is the list of snps obtained from varscan.
    """
    bam_name = find_chromosome_name(bam_file, samtools)
    sequence = list(SeqIO.parse(reference, "fasta"))[0].seq
    length_ref = len(sequence)
    bam_start = 0
    bam_end = length_ref - 1
    shelve_name = bam_file.split(".bam")[0] + ".db"

    if shelve_name not in os.listdir(os.getcwd()):
        database_sam(bam_file, bam_name, bam_start, bam_end, shelve_name)

    snp_list = [int(line.split("\t")[1]) for line in snp_report[1:]]
    if not snp_list:
        return False
    if len(snp_list) == 1:
        return True
    nbr_snps = float(len(snp_list))
    db = shelve.open(shelve_name)
    snp_pairs = [(snp_list[x], snp_list[x + 1]) for x in range(len(snp_list) - 1)]
    snp_pairs_found = count_overlapping_reads(snp_pairs, db)
    db.close()

    if stringency == 0 and snp_pairs_found > 0: #  case 0
        # print("case 0")
        return True
    elif stringency < 1 and (snp_pairs_found) / (nbr_snps - 1) > stringency: #  case 1
        # print("case 1")
        return True
    elif stringency == 1 and snp_pairs_found == len(snp_pairs): #  case 2
        # print("case 2")
        return True
    return False


def fetch_contig_in_bam(record, cfg, temporary_folder, samtools_version):
    """
    For a contig in an alignment, extract its reads 
    and use them to form a bam file that can be used for allele phasing.
    """
    samtools = cfg["samtools_path"]
    bcftools = cfg["bcftools_path"]
    vcfutils = cfg["vcfutils_path"]
    quality = cfg["varscan_min_avg_qual"]
    proper_pairs = cfg["use_only_properly_paired_reads"]
    paired_samples = cfg["transcriptomes_with_paired_reads"]
    read_validity = cfg["check_reads_validity"]
    stringency = cfg["validity_stingency"]
    fasta_list = []

    taxon = record.name.split("|")[0]
    contig = record.name.split("|")[1]
    transcriptome = os.path.join(cfg["bam_folder"], taxon)
    template = "An exception of type {0} occurred when trying to run command {2}. \
    Arguments:\n{1!r}"
    if paired_samples == ["ALL"]:
        paired_samples = [taxon]

    #  Extract the contig in a sam alignment, not a bam file, in order to
    #   to keep the headers of all contigs in the transcriptome.
    cmd = "{samtools} view -o {0}_WITH_{1}.sam {2}_sorted.bam {1}".format(
        taxon, contig, transcriptome, samtools=samtools
    )
    try:
        run_subprocess(cmd)
    except Exception as ex:
        message = template.format(type(ex).__name__, ex.args, cmd)
        print(message)
        sys.exit(message)

    #  Convert sam to bam using an indexed reference taxon.
    cmd = "{samtools} view -bt {0}_WITH_{1}_ref.fas.fai {0}_WITH_{1}.sam \
> {0}_WITH_{1}.bam".format(
        taxon, contig, samtools=samtools
    )

    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
    )
    stdout_str, stderr_str = process.communicate()

    def taxon_in_list(taxon, list_taxa):
        for item in list_taxa:
            if taxon in item:
                return True

    if process.returncode == 0:
        #  Remove improper pairs if requested.
        if proper_pairs and taxon_in_list(taxon, paired_samples):
            prefix = "{0}_WITH_{1}".format(taxon, contig)
            os.rename(prefix + ".bam", prefix + "_temp.bam")
            cmd = "{samtools} view -b -f 0x0002 {prefix}_temp.bam \
> {prefix}.bam".format(
                samtools=samtools, prefix=prefix
            )
            # print("Removing improper pairs.")
            try:
                run_subprocess(cmd)
            except:
                os.rename(prefix + "_temp.bam", prefix + ".bam")

            #  Check for empty bam file.
            cmd = [samtools, "depth", "-q", "30", prefix + ".bam"]
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
            out, err = process.communicate()
            output = out.decode().split("\n")

            if not output:
                fasta_list.append("{0}_WITH_{1}_ref.fas".format(taxon, contig))
                return fasta_list

        #  Sort and index the bam file.
        if samtools_version == 0:
            cmd = "{samtools} sort {0}_WITH_{1}.bam {0}_WITH_{1}_sorted".format(
                taxon, contig, samtools=samtools
            )
        elif samtools_version == 1:
            cmd = "{samtools} sort -T {temp}/{0}_WITH_{1}.bam {0}_WITH_{1}.bam\
 -o {0}_WITH_{1}_sorted.bam".format(
                taxon, contig, samtools=samtools, temp=temporary_folder
            )
        try:
            run_subprocess(cmd)
        except Exception as ex:
            message = template.format(type(ex).__name__, ex.args, cmd)
            print(message)
            sys.exit(message)

        cmd = "{samtools} index {0}_WITH_{1}_sorted.bam".format(
            taxon, contig, samtools=samtools
        )
        try:
            run_subprocess(cmd)
        except Exception as ex:
            message = template.format(type(ex).__name__, ex.args, cmd)
            print(message)
            sys.exit(message)

        #  Call snps.
        call_snp(record, cfg)
        snp_report = []
        #  Tag that determine whether there are some reads for both allelic
        #  variants and not just indels since the phasing does not work
        #  on indels.
        position_validity = False
        for line in open("{0}_WITH_{1}_spns.txt".format(taxon, contig), "rb"):
            snp_report.append(line.decode())
        for line in snp_report[1:]:
            parsed_line = line.split("\t")
            if parsed_line[4] != "0" and parsed_line[5] != "0":
                position_validity = True

        #  Check reads validity if requested (read_validity).
        validity_switch = True
        if read_validity:
            validity_switch = check_validity_stringency(
                samtools,
                taxon + "_WITH_" + contig + "_sorted.bam",
                stringency,
                taxon + "_WITH_" + contig + "_ref.fas",
                snp_report,
            )
            print("validity check is: {}".format(validity_switch))

        #  Phase the allele only is there are some real snps and that
        #  the alignment has passe the validity check.
        if len(snp_report) > 1 and position_validity and validity_switch:
            if samtools_version == 0:
                cmd = "{samtools} phase -A -F -Q {quality} -D 10000 \
{0}_WITH_{1}_sorted.bam -o {0}_WITH_{1}_phased".format(
                    taxon, contig, samtools=samtools, quality=quality
                )
            elif samtools_version == 1:
                cmd = "{samtools} phase -A -F -Q {quality} -D 10000 \
{0}_WITH_{1}_sorted.bam -b {0}_WITH_{1}_phased".format(
                    taxon, contig, samtools=samtools, quality=quality
                )

            try:
                stdout = run_subprocess(cmd, get_stdout=True)
            except Exception as ex:
                message = template.format(type(ex).__name__, ex.args, cmd)
                print(message)
                sys.exit(message)
            
            # Check if samtools was able to successfully phase
            if not stdout or stdout.endswith(b"#errors0  #supp1  #err1\nCC\nCC\n"):
                fasta_list.append("{0}_WITH_{1}_ref.fas".format(taxon, contig))
                return fasta_list
            #  Case not two alleles where produced.
            if "{0}_WITH_{1}_phased.0.bam".format(taxon, contig) not in os.listdir(
                os.getcwd()
            ) or "{0}_WITH_{1}_phased.1.bam".format(taxon, contig) not in os.listdir(
                os.getcwd()
            ):
                sys.exit("[error]: samtools phase has failed to produce two alleles")

            #  Case phasing was succesful.
            for i in [0, 1]:
                if samtools_version == 0:
                    cmd = "{samtools} sort {0}_WITH_{1}_phased.{2}.bam \
{0}_WITH_{1}_phased.{2}_sorted".format(
                        taxon, contig, i, samtools=samtools
                    )
                elif samtools_version == 1:
                    cmd = "{sam} sort -T {temp}/{0}_WITH_{1}_phased.{2}_sorted\
 {0}_WITH_{1}_phased.{2}.bam -o {0}_WITH_{1}_phased.{2}_sorted.bam".format(
                        taxon, contig, i, sam=samtools, temp=temporary_folder
                    )
                try:
                    run_subprocess(cmd)
                except Exception as ex:
                    message = template.format(type(ex).__name__, ex.args, cmd)
                    print(message)
                    sys.exit(message)

                cmd = "{sam} index {0}_WITH_{1}_phased.{2}_sorted.bam".format(
                    taxon, contig, i, sam=samtools
                )
                try:
                    run_subprocess(cmd)
                except Exception as ex:
                    message = template.format(type(ex).__name__, ex.args, cmd)
                    print(message)
                    sys.exit(message)

                #  Get fastq consensuses from phased bam files.
                if samtools_version == 0:
                    cmd = "{sam} mpileup -uf {0}_WITH_{1}_ref.fas \
{0}_WITH_{1}_phased.{2}_sorted.bam | \
{bcf} view -cg  -D -d1 - | {vcf} vcf2fq > {0}_WITH_{1}_phased.{2}.fq".format(
                        taxon, contig, i, sam=samtools, bcf=bcftools, vcf=vcfutils
                    )

                elif samtools_version == 1:
                    cmd = "{sam} mpileup -uf {0}_WITH_{1}_ref.fas \
{0}_WITH_{1}_phased.{2}_sorted.bam | \
{bcf} call -c | {vcf} vcf2fq > {0}_WITH_{1}_phased.{2}.fq".format(
                        taxon, contig, i, sam=samtools, bcf=bcftools, vcf=vcfutils
                    )
                try:
                    run_subprocess(cmd)
                except Exception as ex:
                    message = template.format(type(ex).__name__, ex.args, cmd)
                    print(message)
                    sys.exit(message)

                fastq_product = "{0}_WITH_{1}_phased.{2}.fq".format(taxon, contig, i)

                if (
                    fastq_product not in os.listdir(os.getcwd())
                    or os.stat(fastq_product).st_size < 4
                ):
                    s = "[error]: fastq file {0} has not been generated \n\
                    with command {1}".format(
                        fastq_product, cmd
                    )
                    print(s)
                    sys.exit("line 535, " + s)

                #  Convert fastq to fasta
                if (
                    fastq_product in os.listdir(os.getcwd())
                    and os.stat(fastq_product).st_size > 4
                ):
                    basename = "{0}_WITH_{1}_phased.{2}".format(taxon, contig, i)
                    SeqIO.convert(
                        basename + ".fq", "fastq", basename + ".fasta", "fasta"
                    )
                    fasta_list.append(basename + ".fasta")

                else:
                    if i == 0:
                        # fasta_list.append("{0}_WITH_{1}_ref.fas".format(taxon, contig))
                        sys.exit("line 551")
                        break
            return fasta_list

        else:
            fasta_list.append("{0}_WITH_{1}_ref.fas".format(taxon, contig))
            return fasta_list

    else:
        fasta_list.append("{0}_WITH_{1}_ref.fas".format(taxon, contig))
        return fasta_list


def assemble_final_fasta(final_fasta_list, name_alignment):
    """
    Assemble the phased sequences into complete alignments.
    """
    phased_alignment_name = name_alignment + "_phased.fasta"
    fasta_new_records = []

    for item in final_fasta_list:
        if "phased" in item:
            splitted = item[:-15].split("_WITH_")
            taxon = splitted[0]
            contig = splitted[1]
            i = item[-7]
            rec = list(SeqIO.parse(item, "fasta"))
            new_record = SeqRecord(
                seq=Seq(
                    str(rec[0].seq).replace("n", "").replace("-", ""),
                    IUPAC.unambiguous_dna,
                ),
                name=taxon + "|" + contig + "_allele_{}".format(i),
                id=taxon + "|" + contig + "_allele_{}".format(i),
                description="",
            )
            fasta_new_records.append(new_record)
        elif "_ref.fas" in item:
            splitted = item[:-8].split("_WITH_")
            taxon = splitted[0]
            contig = splitted[1]
            rec = list(SeqIO.parse(item, "fasta"))
            new_record = SeqRecord(
                rec[0].seq,
                name=taxon + "|" + contig,
                id=taxon + "|" + contig,
                description="",
            )
            fasta_new_records.append(new_record)
        else:
            splitted = item.split("_WITH_")
            taxon = splitted[0]
            contig = splitted[1].split(".fa")[0]
            rec = list(SeqIO.parse(item, "fasta"))
            new_record = SeqRecord(
                rec[0].seq,
                name=taxon + "|" + contig,
                id=taxon + "|" + contig,
                description="",
            )
            fasta_new_records.append(new_record)
    SeqIO.write(fasta_new_records, phased_alignment_name, "fasta")


def errors(msg, logger):
    """
    Log error and raise exception.
    """
    logger.error(msg)
    raise Exception(msg)


def check_settings(cfg, logger):
    """
    Assertain that the parameters in the config file are correctly set.
    """
    sanity_checker(cfg, logger)
    varscan_min_fq = cfg["varscan_min_fq"]
    if not 0 <= varscan_min_fq <= 1:
        errors(
            "[error]: Wrong value for varscan_min_fq. \n\
               Check that 0 <= varscan_min_fq <= 1"
        )
    varscan_min_avg_qual = cfg["varscan_min_avg_qual"]
    if not 0 < varscan_min_avg_qual <= 40:
        errors(
            "[error]: Wrong value for varscan_min_avg_qual. \n\
               Check that 0 < varscan_min_avg_qual <= 40"
        )

    if cfg["use_only_properly_paired_reads"]:
        paired_transcriptomes = cfg["transcriptomes_with_paired_reads"]
        if paired_transcriptomes == "ALL":
            L = [
                x.split(".bam")[0]
                for x in os.listdir(cfg["bam_folder"])
                if ".bam" in x and ".bai" not in x
            ]
            cfg["transcriptomes_with_paired_reads"] = L
        elif not paired_transcriptomes:
            errors(
                "[error] Some transcriptomes are paired but \
                   the corresponding list 'transcriptomes_with_paired_reads' is empty."
            )

    if cfg["check_reads_validity"]:
        stingency = cfg["validity_stingency"]
        if not 0 <= stingency <= 1:
            errors(
                "[error]: Wrong value for validity_stingency. \n\
                   Check that 0 <= validity_stingency <= 1"
            )


def phasing_contig(fasta, cfg, **kargs):
    """
    Main function for phasing.
    """
    temp = kargs["temp"]
    samtools_version = kargs["samtools_version"]
    excluded_taxa = cfg["taxa_excluded_from_phasing"]
    os.chdir(cfg["input_folder"])

    records_initial = list(SeqIO.parse(fasta, "fasta"))
    new_temp = fasta.split(".")[0]
    os.chdir(temp)
    if new_temp not in [x for x in os.listdir(os.getcwd())]:
        os.mkdir(new_temp)
    os.chdir(new_temp)
    final_fasta_list = []

    for record in records_initial:
        #  If the record is not part of the excluded taxa proceed
        #   with phasing, otherwise return the unphased sequence
        if record.name.split("|")[0] not in excluded_taxa:
            extract_reference(record, cfg)
            #print("extracting reference for", record.name)
            fasta_list = fetch_contig_in_bam(record, cfg, temp, samtools_version)
            final_fasta_list.extend(fasta_list)

        else:
            #print("record {} in excluded taxa".format(record.name))
            taxon = record.name.split("|")[0]
            contig = record.name.split("|")[1]
            new_seq = Seq(
                str(record.seq).replace("n", "").replace("-", ""),
                IUPAC.unambiguous_dna,
            )
            new_record = SeqRecord(new_seq, name=contig, id=contig, description="")
            SeqIO.write(new_record, taxon + "_WITH_" + contig + ".fas", "fasta")
            final_fasta_list.append(taxon + "_WITH_" + contig + ".fas")

    # print("final_fasta_list", final_fasta_list)
    name_align = fasta.split(cfg["input_suffix"])[0]
    assemble_final_fasta(final_fasta_list, name_align)
    shutil.copy(name_align + "_phased.fasta", cfg["output_folder"])

    os.chdir(os.pardir)
    s_i = "Done phasing fasta: {}".format(fasta)
    s_d = "Phasing fasta: {} was successful".format(fasta)
    return (s_i, s_d)

