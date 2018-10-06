# Takes paired-end STAR output Chimeric.out.junction
# Finds "mate-outside" and "mate-inside", with GT/AG and unique mapped reads, makes bigBed files for UCSC genome browser

# Imports
import subprocess, sys, os, datetime
from collections import defaultdict, OrderedDict

import pandas as pd
import numpy as np
from interval import interval
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ptes.lib.general import init_file, writeln_to_file, shell_call, print_time
from ptes.constants import PTES_logger
from ptes.ptes import get_read_interval, one_interval, get_interval_length, get_subseq
from ptes.ucsc.ucsc import order_interval_list, list_to_dict, get_track_list

### Arguments
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="STAR output, Chimeric.out.junction")                      
parser.add_argument("-o","--output", type=str,
                    help="Output folder for results")  
parser.add_argument("-g","--genome", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',    
                    help="Absolute path to genome file")  
parser.add_argument("-t","--tag", type=str,
                    help="Tag name for grouping results, i.e. ENCODE id")                    
args = parser.parse_args()

# Functions
def split_by_p(read_dict):
    '''
    Takes read_dict (OrderedDict),
    returns list of dicts: before p and after p
    '''
    items = read_dict.items()
    for i, item in enumerate(items):
        if 'p' in item[0]:
            return [OrderedDict(items[:i]), OrderedDict(items[(i+1):])]
    return [read_dict]     
            

# Main

chimeric_file = args.input
input_name = chimeric_file+'.filtered'
path_to_file = args.output.rstrip('/')
with open(input_name, 'r') as input_file:
    for line in input_file:
        line_list = line.strip(None).split('\t')
        chrom = line_list[0]
        intron1 = line_list[1]    #donor splice site coord
        chain = line_list[2]
        intron2 = line_list[4]    #acceptor splice site coord
        junction_type = line_list[6]   #junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC        
        read_name = line_list[9]
        coord1 = int(line_list[10])
        cigar1 = line_list[11]
        coord2 = int(line_list[12])
        cigar2 = line_list[13]
        chim_part1 = get_read_interval(cigar1, coord1)   # not mates, chimeric parts!
        chim_part2 = get_read_interval(cigar2, coord2)
        type = None
        if '-' in cigar1 or '-' in cigar2:    #if p<0 than interval for p does not make any sense
            check_p = False            
        else:
            check_p = True    
        