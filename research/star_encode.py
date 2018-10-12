# Takes paired-end STAR output Chimeric.out.junction
# Finds "mate-outside" and "mate-inside", with GT/AG and unique mapped reads, 
# Counts x / (x + y) where x is N(mate-outside) and y is N(mate-inside)
# Makes bigBed files for UCSC genome browser

# Imports
import os
from collections import defaultdict, OrderedDict

from interval import interval
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ptes.lib.general import init_file, writeln_to_file, write_to_file, shell_call
from ptes.constants import PTES_logger
from ptes.ptes import get_read_interval, one_interval, get_interval_length, get_subseq, split_by_p
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

def dict_to_interval(read_dict):
    '''
    Takes read_dict (OrderedDict),
    returns interval of intervals for all features that consume reference
    '''
    output_interval = interval()
    for item in read_dict.items():
        feature = item[0]
        if 'M' in feature or 'D' in feature or 'N' in feature:
            output_interval = output_interval | item[1]
    return output_interval

# Main

# Exons GTF to junctions dict

PTES_logger.info('Reading GTF...')
gtf_exons_name = '/uge_mnt/home/sunnymouse/Human_ref/hg19_exons_prot_coding.gtf'
gtf_donors = defaultdict(set)
gtf_acceptors = defaultdict(set)
with open(gtf_exons_name, 'r') as gtf_exons_file:
    for line in gtf_exons_file:
        line_list = line.strip().split()
        chrom = line_list[0]
        strt = int(line_list[3])
        end = int(line_list[4])
        chain = line_list[6]
        if chain == '+':
            gtf_donors[chrom].add(end+1)
            gtf_acceptors[chrom].add(strt-1)
        elif chain == '-':
            gtf_donors[chrom].add(strt-1)
            gtf_acceptors[chrom].add(end+1)

PTES_logger.info('Reading GTF... done')  

# Reading filtered STAR output
PTES_logger.info('Reading STAR output...')
input_name = args.input
path_to_file = args.output.rstrip('/')
outside_name = 'mate_outside.junction'
outside_list = []
init_file(outside_name, folder = path_to_file)

mates_inside = 0
mates_outside = 0
annot_donors = 0
annot_acceptors = 0
with open(input_name, 'r') as input_file:
    for line in input_file:
        line_list = line.strip(None).split('\t')
        chrom = line_list[0]
        donor_ss = int(line_list[1])    #donor splice site coord
        chain = line_list[2]
        acceptor_ss = int(line_list[4])    #acceptor splice site coord
        junction_type = line_list[6]   #junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC        
        read_name = line_list[9]
        coord1 = int(line_list[10])
        cigar1 = line_list[11]
        coord2 = int(line_list[12])
        cigar2 = line_list[13] 
        chim_part1 = get_read_interval(cigar1, coord1)   # not mates, chimeric parts!
        chim_part2 = get_read_interval(cigar2, coord2)
        if 'p' in cigar1:
            splits = split_by_p(chim_part1)
            if chain == '+':
                mate1 = dict_to_interval(splits[0])
                mate_intervals = dict_to_interval(splits[1]) | dict_to_interval(chim_part2)
                mate2 = one_interval(mate_intervals)
            if chain == '-':
                mate_intervals = dict_to_interval(splits[0]) | dict_to_interval(chim_part2)
                mate1 = one_interval(mate_intervals)
                mate2 = dict_to_interval(splits[1])                
        elif 'p' in cigar2:
            splits = split_by_p(chim_part2)
            if chain == '+':
                mate_intervals = dict_to_interval(chim_part1) | dict_to_interval(splits[0])
                mate1 = one_interval(mate_intervals)
                mate2 = dict_to_interval(splits[1])
            if chain == '-':   
                mate_intervals = dict_to_interval(chim_part1) | dict_to_interval(splits[1])
                mate1 = one_interval(mate_intervals)
                mate2 = dict_to_interval(splits[0])
        print
        print cigar1, cigar2
        print coord1, coord2
        print chain
        print chim_part1
        print chim_part2
        print
        print mate_intervals
        print mate1
        print mate2
        mate_intersection = mate1 & mate2
        if mate_intersection == interval():   # zero intersection
            print 'mate outside'
            mates_outside += 1  
            outside_list.append(line)
        else:    
            print 'mate inside'
            print 'Length of intersection: %i' % get_interval_length(mate_intersection)
            mates_inside += 1        
        if donor_ss in gtf_donors:
            annot_donors += 1
        if acceptor_ss in gtf_acceptors:            
            annot_acceptors += 1

 
PTES_logger.info('Reading STAR output... done')  
print 'Inside: %i' % mates_inside      
print 'Outside: %i' % mates_outside      
print 'Annot donors: %i' % annot_donors      
print 'Annot acceptors: %i' % annot_acceptors   

writeln_to_file(''.join(outside_list), outside_name, folder = path_to_file)   