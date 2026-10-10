#!/usr/bin/env python3
# _*_ coding: UTF-8 _*_

#command
#python ./count-stop_v86.py -i ${file) -r /reference_geneme/genome.fa --StrandSelect fwd > ${file}_stop
#
import pysam
import argparse
from Bio import SeqIO
import sys

# Version information START --------------------------------------------------
VERSION_INFO = \
"""
    Author: MENG Howard

    Version-01:
      2018-03-09, set up the code

    Version-02:
      2018-03-13, fix the code that can run with all site,
      include the reads without CMC stop

    Version-03:
      2018-03-15, random access the Ref FASTA;

    Version-04
      2018-03-28, use Bio package to access the Ref FASTA;

    Version-05
      2018-03-28, add load method select; use memory efficiency method to run the code

    Version-06
      2018-04-25, fix the reverse strand analysis bug!!!
          1. only support the [fwd] or [rev] mode
          2. the reverse code have a huage bug, but I fixed it.

    E-Mail: meng_howard@126.com
"""
# Version information END ----------------------------------------------------

# Learning Part START --------------------------------------------------------
LEARNING_PART = \
"""
    1.pysam alignment result is 0-base pythonic index, for example a reads with
        alignment index (4,25,21M), means the 5th base of reference is the start
        base but not contain the 26th base of reference.
"""

# Learning Part END-----------------------------------------------------------

###############################################################################
# def function
###############################################################################
def make_fa_dict(fasta_file_path):
    """
    <fasta_file_path>
        fasta file path
    <return>
        A fasta dict, key = fa name, value = fa seq
    """
    fasta_dict = {}
    for seq_record in SeqIO.parse(fasta_file_path,"fasta"):
        fasta_dict_get = fasta_dict.get(seq_record.id)
        if fasta_dict_get is None:
            fasta_dict[seq_record.id] = str(seq_record.seq)
        else:
            if len(fasta_dict_get) < len(str(seq_record.seq)):
                fasta_dict[seq_record.id] = str(seq_record.seq)
    return(fasta_dict)


def reverse_complement_seq(seq="",reverse_state=True,complement_state=True):
    """
    INPUT:
        <seq>
            A DNA sequence like 'ATCG'
        <reverse_state>
            If return a reverse sequence, default if [True]
        <complement_state>
            If return a sequence with complementary information [True]
    RETURN:
        A sequence.
    """
    if reverse_state:
        seq = seq[::-1]

    if complement_state:
        reverse_complement_table = str.maketrans("ATCG", "TAGC")
        seq = seq.translate(reverse_complement_table)

    return(seq)


def echo_site_info(ref_name,ref_seq,current_position,stop_count_list,coverage_count_list,strand_select = "fwd"):
    """
    INPUT:
        <ref_name>
            chromesome name like "chr1"
        <ref_seq>
            chromesome sequence as string format
        <current_position>
            1-based current_position, means the start index of site_info
        <stop_count_list>
            stop count list, length = BLOCK_SIZE
        <coverage_count_list>
            coverage count list length = BLOCK_SIZE

    RETURN:
        <None>

    HELP:
        This function will return None object, and print the site information on the stdout.
    """
    # output format
    output_fmt = "{ref_name},{stop_position},{strand},{base_info},{site_index},{stop_reads},{total_reads},{stop_rate}"

    # strand info
    strand_info = "+"
    if strand_select == "rev":
        strand_info = "-"

    # get chrom_index
    if current_position == 0:
        chrom_start_index = 0
        chrom_end_index =  len(stop_count_list)
    else:
        chrom_start_index = current_position - 1
        chrom_end_index = current_position + len(stop_count_list) - 1

    # get real base
    chrom_seq_fwd = ref_seq[chrom_start_index:chrom_end_index].upper()
    chrom_seq_real_base = chrom_seq_fwd
    if strand_info == "-":
        chrom_seq_real_base = reverse_complement_seq(chrom_seq_fwd,reverse_state=False,complement_state=True)

    # output site information
    for index,cover_count in enumerate(coverage_count_list):
        # pass the first base
        if (index == 0) and (current_position == 0):
            continue

        total_count = cover_count + stop_count_list[index]

        if (total_count > 0) and (index+1 <= len(chrom_seq_fwd)):
            # make site index(1-base)
            site_index = "%s_%d_%s_%s" % (ref_name,current_position + index,strand_info,chrom_seq_fwd[index])

            # calculate stop rate
            stop_rate = round(stop_count_list[index] / 1.0 / total_count * 100,2)

            # make output str
            site_str = output_fmt.format(ref_name = ref_name,
                                       stop_position = current_position + index,
                                       strand = strand_info,
                                       base_info = chrom_seq_real_base[index],
                                       site_index = site_index,
                                       stop_reads = stop_count_list[index],
                                       total_reads = total_count,
                                       stop_rate = stop_rate)
            print(site_str)

    return(None)

###############################################################################
# read parameters
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="This tools can count stop information from BAM file.")

    parser.add_argument("-i", "--InputBam",
        help="BAM file, which have to be sorted",required=True)

    parser.add_argument("-r", "--Ref",
        help="Reference with FASTA format",required=True)

    parser.add_argument("-c", "--Coverage",
        help="Site coverage threshold, default = 1",default="1")

    parser.add_argument("--StrandSelect",
        help="Coutn site information use both starnd information,\
        [fwd] means select + ,\
        [rev] means select -, default = [fwd]",
        default = "fwd")

    parser.add_argument("--BlockSize",
        help="The block size of counting job, the larger number means the larger memory space costs,\
        default=[1000000]",default="1000000")

    parser.add_argument("--LargeGenome",
        help="If set [Yes] will take more running time but low memory requirment. default=[No]",default="No")

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # load the paramters
    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ARGS = parser.parse_args()

    # coverage threshold
    OUTPUT_COVER_CUTOFF = eval(ARGS.Coverage)

    # get select information
    STRAND_SELECT = ARGS.StrandSelect
    if (STRAND_SELECT != "fwd") and (STRAND_SELECT != "rev"):
        raise IOError("--StrandSelect have to set one of [fwd] or [rev].")

    # input BAM file
    INPUT_FILE = pysam.AlignmentFile(ARGS.InputBam)

    # block size
    BLOCK_SIZE = int(ARGS.BlockSize)

    # use this method to access the fasta file
    if ARGS.LargeGenome == "Yes":
        REF_SEQ_DICT = SeqIO.index(ARGS.Ref,"fasta")
    else:
        REF_SEQ_DICT = make_fa_dict(ARGS.Ref)

    ##########################################################################
    # initialize the paramters
    ##########################################################################
    # all code write as 1-based system !
    # ref parameter
    current_ref_name = ""
    current_ref_seq = ""
    current_ref_len = 0

    # record some factors
    current_position = 0
    cur_coverage_count_list = []
    cur_stop_count_list = []

    for read in INPUT_FILE:
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # 1.check reads basic state information
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # get the strand information
        if read.is_reverse:
            strand_info = "-"
        else:
            strand_info = "+"

        # the reference_index is error
        if read.reference_id < 0:
            continue

        # the reference_name not match the ref_seq_dict
        if not read.reference_name in REF_SEQ_DICT.keys():
            continue

        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # 2.fwd and rev with different strategy
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        # select reads to process
        if (STRAND_SELECT == "fwd" and strand_info == "+") or (STRAND_SELECT == "rev" and strand_info == "-"):
            if STRAND_SELECT == "fwd":
                # 2.1 initialize, set the paramters >>>>>>>>>>>>>>>>>>>>>>>>>>>
                if current_ref_name == "":
                    # initialize: load FASTA sequence
                    current_ref_name = read.reference_name

                    if ARGS.LargeGenome == "Yes":
                        current_ref_seq = str(REF_SEQ_DICT.get(current_ref_name).seq)
                    else:
                        current_ref_seq = REF_SEQ_DICT.get(current_ref_name)

                    # initialize: make list
                    current_list_len = BLOCK_SIZE
                    current_position = read.reference_start

                    cur_coverage_count_list = [0] * current_list_len
                    cur_stop_count_list = [0] * current_list_len

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # 2.2 change the ref_name
                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                elif current_ref_name != read.reference_name:
                    # get output
                    echo_site_info(ref_name=current_ref_name,
                                ref_seq=current_ref_seq,
                                current_position=current_position,
                                stop_count_list=cur_stop_count_list,
                                coverage_count_list=cur_coverage_count_list,
                                strand_select = "fwd")

                    # re-initialize
                    current_ref_name = read.reference_name

                    if ARGS.LargeGenome == "Yes":
                        current_ref_seq = str(REF_SEQ_DICT.get(current_ref_name).seq)
                    else:
                        current_ref_seq = REF_SEQ_DICT.get(current_ref_name)

                    # initialize: make list
                    current_list_len = BLOCK_SIZE
                    current_position = read.reference_start

                    cur_coverage_count_list = [0] * current_list_len
                    cur_stop_count_list = [0] * current_list_len

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # 2.3 out the block region
                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                elif read.reference_end > (current_position + BLOCK_SIZE -1):
                    # get output
                    echo_site_info(ref_name=current_ref_name,
                                ref_seq=current_ref_seq,
                                current_position=current_position,
                                stop_count_list=cur_stop_count_list,
                                coverage_count_list=cur_coverage_count_list,
                                strand_select = "fwd")

                    # re-initialize: make list
                    current_list_len = BLOCK_SIZE
                    current_position = read.reference_start

                    cur_coverage_count_list = [0] * current_list_len
                    cur_stop_count_list = [0] * current_list_len

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # 2.4 count the coverage AND stop_info
                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # NOTE (local fix): a read may span a huge intron (> BLOCK_SIZE),
                # so positions beyond the current block window are SKIPPED
                # instead of splitting the read across blocks. This avoids
                # IndexError for position_index - current_position + 1 out of
                # range. Positions inside the window are counted as before.
                for position_index in read.positions:
                    pos_idx = position_index - current_position + 1
                    if 1 <= pos_idx < len(cur_coverage_count_list):
                        cur_coverage_count_list[pos_idx] += 1

                # count stop information
                if 0 <= (read.reference_start - current_position) < len(cur_stop_count_list):
                    cur_stop_count_list[read.reference_start - current_position ] += 1


            elif STRAND_SELECT == "rev":
                # 3.1 initialize, set the paramters >>>>>>>>>>>>>>>>>>>>>>>>>>>
                if current_ref_name == "":
                    # initialize: load FASTA sequence
                    current_ref_name = read.reference_name

                    if ARGS.LargeGenome == "Yes":
                        current_ref_seq = str(REF_SEQ_DICT.get(current_ref_name).seq)
                    else:
                        current_ref_seq = REF_SEQ_DICT.get(current_ref_name)

                    # initialize: make list
                    current_list_len = BLOCK_SIZE
                    current_position = read.reference_start + 1

                    cur_coverage_count_list = [0] * current_list_len
                    cur_stop_count_list = [0] * current_list_len

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # 3.2 change the ref_name
                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                elif current_ref_name != read.reference_name:
                    # get output
                    echo_site_info(ref_name=current_ref_name,
                                ref_seq=current_ref_seq,
                                current_position=current_position,
                                stop_count_list=cur_stop_count_list,
                                coverage_count_list=cur_coverage_count_list,
                                strand_select = "rev")

                    # re-initialize
                    current_ref_name = read.reference_name

                    if ARGS.LargeGenome == "Yes":
                        current_ref_seq = str(REF_SEQ_DICT.get(current_ref_name).seq)
                    else:
                        current_ref_seq = REF_SEQ_DICT.get(current_ref_name)

                    # initialize: make list
                    current_list_len = BLOCK_SIZE
                    current_position = read.reference_start + 1

                    cur_coverage_count_list = [0] * current_list_len
                    cur_stop_count_list = [0] * current_list_len

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # 3.3 out the block region
                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                elif read.reference_end >= (current_position + BLOCK_SIZE -1):
                    # get output
                    echo_site_info(ref_name=current_ref_name,
                                ref_seq=current_ref_seq,
                                current_position=current_position,
                                stop_count_list=cur_stop_count_list,
                                coverage_count_list=cur_coverage_count_list,
                                strand_select = "rev")

                    # re-initialize: make list
                    current_list_len = BLOCK_SIZE
                    current_position = read.reference_start + 1

                    cur_coverage_count_list = [0] * current_list_len
                    cur_stop_count_list = [0] * current_list_len

                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # 3.4 count the coverage AND stop_info
                # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                # NOTE (local fix): same as 2.4, positions beyond the current
                # block window (huge-intron reads) are SKIPPED to avoid IndexError.
                for position_index in read.positions:
                    pos_idx = position_index - current_position + 1
                    if 1 <= pos_idx < len(cur_coverage_count_list):
                        cur_coverage_count_list[pos_idx] += 1

                # count stop information
                stop_idx = read.reference_end - current_position + 1
                if 0 <= stop_idx < len(cur_stop_count_list):
                    cur_stop_count_list[stop_idx] += 1
                else:
                    sys.stderr.write("\n" + str(read) + "\n")
                    sys.stderr.write(str(read.reference_end) + "\n")
                    sys.stderr.write(str(current_position) + "\n")
                    sys.stderr.write(str(len(cur_stop_count_list))+ "\n")
                    sys.stderr.write("-"*80 + "\n")
        else:
            continue

    # output the remaining
    if STRAND_SELECT == "fwd" :
        echo_site_info(ref_name=current_ref_name,
                    ref_seq=current_ref_seq,
                    current_position=current_position,
                    stop_count_list=cur_stop_count_list,
                    coverage_count_list=cur_coverage_count_list,
                    strand_select = "fwd")

    elif STRAND_SELECT == "rev":
        echo_site_info(ref_name=current_ref_name,
                    ref_seq=current_ref_seq,
                    current_position=current_position,
                    stop_count_list=cur_stop_count_list,
                    coverage_count_list=cur_coverage_count_list,
                    strand_select = "rev")























#
