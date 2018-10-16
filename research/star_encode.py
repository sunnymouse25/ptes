# Takes paired-end STAR output Chimeric.out.junction
# Finds "mate-outside" and "mate-inside", with GT/AG and unique mapped reads, 
# Counts x / (x + y) where x is N(mate-outside) and y is N(mate-inside)
# Copies lines with mates outside to mates_outside.junction to make bigBed files for UCSC genome browser

# Imports

### Arguments
import argparse

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file
from ptes.ptes import annot_junctions, \
    mate_intersection, return_mates

# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
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


# Exons GTF to junctions dict

PTES_logger.info('Reading GTF...')
gtf_exons_name = '/uge_mnt/home/sunnymouse/Human_ref/hg19_exons.gtf'
gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)

PTES_logger.info('Reading GTF... done')  

# Reading filtered STAR output
PTES_logger.info('Reading STAR output...')
input_name = args.input
path_to_file = args.output.rstrip('/')
outside_name = 'mate_outside.junction'
outside_list = []
init_file(outside_name, folder = path_to_file)

mates = {'inside': 0, 'outside': 0, 'non-chim' : 0}

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
        print 'Length of chim junction %i' % (donor_ss-acceptor_ss)
        if chain == '+':
            if donor_ss < acceptor_ss:
                print 'Not chimeric? ', donor_ss, acceptor_ss
                mates['non-chim'] += 1
                continue
        elif chain == '-':
            if donor_ss > acceptor_ss:
                print 'Not chimeric? ', donor_ss, acceptor_ss
                mates['non-chim'] += 1
                continue
        mate1, mate2 = return_mates(cigar1=cigar1,
                                    coord1=coord1,
                                    cigar2=cigar2,
                                    coord2=coord2,
                                    chain=chain)
        interval_intersection = mate_intersection(mate1, mate2)
        mates[interval_intersection] += 1
        if interval_intersection == 'outside':
            outside_list.append(line)
        if donor_ss in gtf_donors:
            annot_donors += 1
        if acceptor_ss in gtf_acceptors:            
            annot_acceptors += 1

PTES_logger.info('Reading STAR output... done')  
print 'Inside: %i' % mates['inside']
print 'Outside: %i' % mates['outside']
print 'Intron too large: %i' % mates['non-chim']
print 'Annot donors: %i' % annot_donors
print 'Annot acceptors: %i' % annot_acceptors   

writeln_to_file(''.join(outside_list), outside_name, folder = path_to_file)   