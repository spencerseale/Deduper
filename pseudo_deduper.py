#!/usr/bin/env python

the problem:

Due to amplification, a common component of NGS library preparation protocols, duplicate
DNA fragments are produced in the final library. During sequencing, these duplicates
will cause a bias in reads of these duplicates that is a product of amplification known
as PCR bias. In RNA-seq especially, this is a problem when the intention is to compare
expression of samples and genes as the bias can unnaturally skew this data. Therefore,
it is necessary to remove this bias after sequencing and prior to analysis

pcr duplicates contain, in order to be checked by script:

same chromosome
same UMI
same position
same strand

the general strategy:

-open the sorted sam file for reading, open an output file for writing
-step through the input sam file, line by line and compare each current line to the dup_tracker dictionary that tracks position, strand direction, and umi
-run current_line function to extract position, umi, and strand directionality (boolean) from current line
    -in this order:
        -check whether UMI is one of the UMIs used in the experiment
            -if no match, then skip line and move to the next line, this UMI is not real and should not be trusted
        -check if current mapped read is on same chromosome as previous read
            -if new chromsome is reached compared to previous line, wipe dup_tracker and proceed
        -run pos_adjuster function to adjust the position of the current read
            -takes in account I and D if reverse strand is being analyzed
                -if I in cigar of a reverse read, subtract number before the I from pos in this_line tuple
                -if D in cigar of a reverse read, add number before the D to pos in this_line tuple
            -adjusts for soft clipping from end of cigar if reverse strand, looks only at beginning of cigar for soft clipping if forward read
                -for forward, if S at the start of cigar, subtract number before S from pos in this_line tuple
                -for reverse, if S at the end of cigar, subtract number before S from pos in this_line tuple
        -Check this_line tuple if in dictionary, if not in dictionary, add to dup_tracker and write out line
            -if in dictionary, clear this_line tuple and move on to next line in the sam file.

#constants

#load umi into a dictionary, this is necessary to check for the umi's in the sam file
umi_dict = {}
with open("./STL96.txt", "r") as umiopen:
    for line in umiopen:
        umi_dict[line.strip()] = line.strip()

#functions

#function to take in user input
def get_args():
    """
    take user input, no examples necessary
    """
    parser = argparse.ArgumentParser(description="remove PCR duplicates from fastq files")
    parser.add_argument("-d", "--dir", help="Input directory containing sam file to be read", type=str)
    parser.add_argument("-s", "--sam", help="Input sam file to be analyzed for PCR duplicates", type=str)
    return parser.parse_args()

def write_read(line, pointer):
    """
    write current read out to specified file using file pointer and current line
    example input: K00337:113:HN7KGBBXX:4:1116:22983:39629:AACGCCAT	99	10	33369999	255	90M	=	33370333	530	GCCAGCAATGACACAAGTCCTTTTTTGTCAAGTTCAATGAGCTTAAGGTGCTGTAACGTTTTCAGTGCAGAATCGTTGACTTTCATGTTT	JJJJJJJJJJJJJJJFJJFFJJJJJJJJFJFFJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFA<<FJJJJJJJJJJJJJJJJ	NH:i:1	HI:i:1	AS:i:182	nM:i:0
    example output: K00337:113:HN7KGBBXX:4:1116:22983:39629:AACGCCAT	99	10	33369999	255	90M	=	33370333	530	GCCAGCAATGACACAAGTCCTTTTTTGTCAAGTTCAATGAGCTTAAGGTGCTGTAACGTTTTCAGTGCAGAATCGTTGACTTTCATGTTT	JJJJJJJJJJJJJJJFJJFFJJJJJJJJFJFFJJJFJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJFA<<FJJJJJJJJJJJJJJJJ	NH:i:1	HI:i:1	AS:i:182	nM:i:0
    """
    pointer.write(line, "w")
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
    parts = line.strip().split("\t")
    rev = False
    #clecking bitwise flag in column 2
    if ((parts[1] & 16) = 16):
        rev = True
    this_line = (parts[0][len(parts[0])-8:], parts[3], rev)
    return this_line

#not completely sure how this function will accomplish goal
#will check directionality
#if reverse, ajust for I and D and look for S at the end of cigar string
#if forward, will only look for S at the start of cigar string
def pos_adjuster(this_line):
    """
    adjusts index 1 of this_line tuple to true position by accounting for muliple factors
    example input: (100, "AATTCCGG", True)
    example output: (98, "AATTCCGG", True)
    """
    if rev:
        check for S at end of cigar
        check for I and check for D
    else:
        check for S at start of cigar
    return this_line






#PSEUDOCODE ENDS HERE






args = get_args()

#sam file will be sorted via sort_sam.srun
directory = args.dir
input_sam = args.sam
output_sam = input_sam[:len(input_sam)-4]+"_deduped.sam"

#setting string to hold previous chromosome
prev_chromo = ""
with open(output_sam, "w") as outsam:
    with open(directory+input_sam, "r") as insam:
        for line in insam:
            #skip lines beginning with @ at the start of sam file
            if not line.startswith("@"):
                this_line = current_line(line)
                this_line = pos_adjuster(this_line)
                #checking if umi in read is one of the umis used in the experiment
                if this_line[1] in umi_dict:
                    #checking chromosome of current read and whether this read has been seen before and is pcr duplicate
                    if line.split()[2] == prev_chromo
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

