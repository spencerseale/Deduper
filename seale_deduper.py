#!/usr/bin/env python

import argparse
import re
# import sys

#functions

#function to take in user input
def get_args():
    """
    take user input, no examples necessary
    """
    parser = argparse.ArgumentParser(description="remove PCR duplicates from fastq files")
    parser.add_argument("-f", "--file", help="Input sam file to be analyzed for PCR duplicates", type=str, required=False)
    parser.add_argument("-u", "--umi", help="File containing UMI sequences used in experiment", type=str, required=False)
    parser.add_argument("-p", "--paired", help="Input sam file is from a paired-end experiment", type=str, required=False)
    #parser.add_argument("-h", "--help", help="Flag to give help message", type=str, required = False)
    return parser.parse_args()

def write_read(line, pointer):
    """
    write current read out to specified file using file pointer and current line
    example input: K00337:113:HN7KGBBXX:4:1116:22983:39629:AACGCCAT	99	10	33369999	255	90M	=	33370333	530	GCCAGCAATGACACAAGTCCTTTTTTGTCAAGTTCAATGAGCTTAAGGTGCTGTAACGTTTTCAGTGCAGAATCGTTGACTTTCATGTTT	JJJJJJJJJJJJJJJFJJFFJJJJJJJJFJFFJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFA<<FJJJJJJJJJJJJJJJJ	NH:i:1	HI:i:1	AS:i:182	nM:i:0
    example output: K00337:113:HN7KGBBXX:4:1116:22983:39629:AACGCCAT	99	10	33369999	255	90M	=	33370333	530	GCCAGCAATGACACAAGTCCTTTTTTGTCAAGTTCAATGAGCTTAAGGTGCTGTAACGTTTTCAGTGCAGAATCGTTGACTTTCATGTTT	JJJJJJJJJJJJJJJFJJFFJJJJJJJJFJFFJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFA<<FJJJJJJJJJJJJJJJJ	NH:i:1	HI:i:1	AS:i:182	nM:i:0
    """
    pointer.write(line)
    return

#this function will write current line to a tuple
#the intention of this function is to track the important contents of current line
def current_line(line):
    """
    stores current line as a 3 component tuple
    checks bitwise flag to see if read being reverse complimented and stores as a boolean
    in order, this list stores UMI, position, rev boolean
    example input: K00337:113:HN7KGBBXX:4:1116:22983:39629:AACGCCAT	99	10	33369999	255	90M	=	33370333	530	GCCAGCAATGACACAAGTCCTTTTTTGTCAAGTTCAATGAGCTTAAGGTGCTGTAACGTTTTCAGTGCAGAATCGTTGACTTTCATGTTT	JJJJJJJJJJJJJJJFJJFFJJJJJJJJFJFFJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFA<<FJJJJJJJJJJJJJJJJ	NH:i:1	HI:i:1	AS:i:182	nM:i:0
    example output: ("AACGCCAT", 33369999, True)
    """
    rev = False
    #clecking bitwise flag in column 2
    if ((int(parts[1]) & 16) == 16):
        rev = True
    this_line = (parts[0][len(parts[0])-8:], parts[3], rev)
    return this_line

#if reverse, move position to where the 5' end is by adding M,N,S, and D
#if forward, will only look for S at the start of cigar string
def pos_adjuster(this_line):
    """
    adjusts index 1 of this_line tuple to true position by accounting for muliple factors
    example input: ("AATTCCGG", 100, True)
    example output: ("AATTCCGG", 98, True)
    """
    this_line = list(this_line)
    cigar = parts[5]
    #if reverse read, need to modify position to be 5', adjust D, M, N, S (only S at end of cigar string)
    if this_line[2]:
        if "M" in cigar:
            M = sum(int(x) for x in (re.findall(r'([0-9]+)M', parts[5])))
            #print(f"M: {M}")
        else:
            M = 0
        if "N" in cigar:
            N = sum(int(x) for x in (re.findall(r'([0-9]+)N', parts[5])))
        else:
            N = 0
        if "D" in cigar:
            D = sum(int(x) for x in (re.findall(r'([0-9]+)D', parts[5])))
        else:
            D = 0
        #reversing order of cigar string to see if S
        if re.match(r'S', cigar[::-1]):
            S = re.findall(r'^S([0-9]+)', cigar[::-1])
            #reversing back to normal order
            S = int(S[0][::-1])
            #print(S)
        else:
            S = 0
        this_line[1] = str(int(this_line[1]) + M + N + D + S)
    #if forward read, only need to check S at start of cigar string to adjust pos
    else:
        if re.match(r'[0-9]+S', cigar):
            S = re.findall(r'^([0-9]+)S', cigar)
            this_line[1] = str(int(this_line[1]) - int(S[0]))
    this_line = tuple(this_line)
    return this_line

args = get_args()

#help statement
# if args.help:
#     print("Due to amplification, a common component of NGS library preparation protocols, duplicate\
#     DNA fragments are produced in the final library. During sequencing, these duplicates\
#     will cause a bias in reads of these duplicates that is a product of amplification known\
#     as PCR bias. In RNA-seq especially, this is a problem when the intention is to compare\
#     expression of samples and genes as the bias can unnaturally skew this data. Therefore,\
#     it is necessary to remove this bias after sequencing and prior to analysis.\
#     \n Currently, this file can only accept single-end sequencing reads.")
#     sys.exit("Error message")


#sam file will be sorted via sort_sam.srun
umi = args.umi
input_sam = args.file
output_sam = input_sam[:len(input_sam)-4]+"_deduped.sam"

#practice files
# input_sam = "./test.sam"
# output_sam = "test_input_deduped.sam"

#load umi into a dictionary, this is necessary to check for the umi's in the sam file
umi_dict = {}
with open(umi, "r") as umiopen:
    for line in umiopen:
        umi_dict[line.strip()] = line.strip()

#setting string to hold previous chromosome
prev_chromo = ""
with open(output_sam, "w") as outsam:
    with open(input_sam, "r") as insam:
        for line in insam:
            #skip lines beginning with @ at the start of sam file
            if not line.startswith("@"):
                parts = line.strip().split("\t")
                this_line = current_line(line)
                this_line = pos_adjuster(this_line)
                #print(this_line)
                #checking if umi in read is one of the umis used in the experiment
                if this_line[0] in umi_dict:
                    #checking chromosome of current read and whether this read has been seen before and is pcr duplicate
                    if line.split()[2] == prev_chromo:
                        if this_line not in dup_tracker:
                            write_read(line, outsam)
                            #add to dictionary
                            dup_tracker[this_line] = this_line
                    else:
                        prev_chromo = line.split()[2]
                        #creates dictionary for first time being seen
                        #clear dictionary everytime a new chromosome is reached
                        dup_tracker = {}
                        write_read(line, outsam)
                        #add to dictionary
                        dup_tracker[this_line] = this_line
            else:
                #write out header lines of sam file
                write_read(line, outsam)

print("--deduping complete--")
