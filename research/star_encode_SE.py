# Takes paired-end STAR output Chimeric.out.junction
# Finds "mate-outside" and "mate-inside", with GT/AG and unique mapped reads,
# Counts x / (x + y) where x is N(mate-outside) and y is N(mate-inside)
# Copies lines with mates outside to mates_outside.junction to make bigBed files for UCSC genome browser

# Imports
from collections import defaultdict


### Arguments
import argparse

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file
from ptes.ptes import annot_junctions, \
    mate_intersection, get_read_interval, dict_to_interval, one_interval, interval_to_string

# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="STAR output, Chimeric.out.junction")
parser.add_argument("-s","--sam", type=str,
                    help="Filtered STAR SAM output, with read_names same as in Chimeric.out.junction")
parser.add_argument("-o","--output", type=str,
                    help="Output folder for results")
parser.add_argument("-g","--genome", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',
                    help="Absolute path to genome file")
parser.add_argument("-t","--tag", type=str,
                    default = 'ENCODE',
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions


# Exons GTF to junctions dict

PTES_logger.info('Reading GTF...')
gtf_exons_name = '/uge_mnt/home/sunnymouse/Human_ref/hg19_exons.gtf'
gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)

PTES_logger.info('Reading GTF... done')

# Reading filtered STAR non-chim output
PTES_logger.info('Reading STAR non-chimeric output...')
sam_name = args.sam
sam_dict = defaultdict(list)
with open(sam_name, 'r') as input_file:
    for line in input_file:
        row = line.strip(None).split('\t')
        read_name = row[0]
        flag = int(row[1])
        chrom = row[2]
        leftpos = row[3]
        cigar = row[5]
        if flag & 16 == 0:
            chain = '+'
        else:
            chain = '-'
        sam_attrs = {
            'chain': chain,
            'chrom' : chrom,
            'leftpos' : leftpos,
            'cigar' : cigar,
        }
        sam_dict[read_name].append(sam_attrs)

PTES_logger.info('Reading STAR non-chimeric output... done')

# Reading filtered STAR output
PTES_logger.info('Reading STAR chimeric output...')
input_name = args.input
path_to_file = args.output.rstrip('/')
outside_name = 'mate_outside.junction'
outside_list = []
outside_intervals = 'mate_outside.test'
init_file(outside_name, folder = path_to_file)
init_file(outside_intervals, folder = path_to_file)
annot_donors = 0
annot_acceptors = 0
mates = {'inside': 0, 'outside': 0, 'non-chim' : 0}
read_names_set = set()

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
        read_names_set.add(read_name)
        if chain == '+':
            if donor_ss < acceptor_ss or abs(donor_ss - acceptor_ss) > 1000000:
                mates['non-chim'] += 1
                continue
        elif chain == '-':
            if donor_ss > acceptor_ss or abs(donor_ss - acceptor_ss) > 1000000:
                mates['non-chim'] += 1
                continue
        chim_part1 = get_read_interval(cigar=cigar1, leftpos=coord1)
        chim_part2 = get_read_interval(cigar=cigar2, leftpos=coord2)
        mate1 = one_interval(dict_to_interval(chim_part1) | dict_to_interval(chim_part2))
        mate2 = None
        for mate in sam_dict[read_name]:
            if mate['cigar'] == cigar1 or mate['cigar'] == cigar2:
                continue
            elif mate['chrom'] == chrom and mate['chain'] != chain:
                mate2 = one_interval(dict_to_interval(get_read_interval(cigar=mate['cigar'], leftpos=mate['leftpos'])))
                interval_intersection = mate_intersection(mate1, mate2)
                mates[interval_intersection] += 1
                if interval_intersection == 'outside':
                    outside_list.append(line)
                    with open(outside_intervals, 'a') as test_file:
                        test_file.write(read_name + '\n')
                        test_file.write(line)
                        test_file.write(interval_to_string(i=mate1)+ '\n')
                        test_file.write(interval_to_string(i=mate2)+ '\n')
                        test_file.write(' '.join(mate.items())+ '\n')
                        test_file.write('\n')
                if donor_ss in gtf_donors:
                    annot_donors += 1
                if acceptor_ss in gtf_acceptors:
                    annot_acceptors += 1

PTES_logger.info('Reading STAR chimeric output... done')

print 'Inside: %i' % mates['inside']
print 'Outside: %i' % mates['outside']
print 'Intron too large: %i' % mates['non-chim']
print 'Annot donors: %i' % annot_donors
print 'Annot acceptors: %i' % annot_acceptors

writeln_to_file(''.join(outside_list), outside_name, folder = path_to_file)

