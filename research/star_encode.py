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
    mate_intersection, return_mates, star_line_dict

# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord
parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="STAR output, Chimeric.out.junction")                      
parser.add_argument("-o","--output", type=str,
                    help="Output folder for results")  
parser.add_argument("-g","--genome", type=str,
                    default = '/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',
                    help="Absolute path to genome file")
parser.add_argument("-gtf","--gtf_annot", type=str,
                    default = '/home/sunnymouse/Human_ref/hg19_exons.gtf',
                    help="Absolute path to genome file")
parser.add_argument("-t","--tag", type=str,
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions


# Exons GTF to junctions dict

PTES_logger.info('Reading GTF...')
gtf_exons_name = args.gtf_annot
gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)

PTES_logger.info('Reading GTF... done')  

# Reading filtered STAR output
PTES_logger.info('Reading STAR output...')
path_to_file = args.output.rstrip('/')
outside_name = 'mate_outside.junction'
outside_list = []
init_file(outside_name, folder = path_to_file)

mates_gtag = {'inside': 0, 'outside': 0, 'non-chim' : 0}
mates_nc = {'inside': 0, 'outside': 0, 'non-chim' : 0}

annot_donors = 0
annot_acceptors = 0
with open(args.input, 'rb') as input_file:
    for line in input_file:
        line_dict = star_line_dict(line=line)
        if line_dict['chrom1'] == line_dict['chrom2'] \
                and line_dict['chain1'] == line_dict['chain2']:
            chrom = line_dict['chrom1']
            chain = line_dict['chain1']
        else:
            PTES_logger.error('Non-filtered STAR output')
            PTES_logger.error('Use awk "$1 ==$4 && $3 ==$6" to filter')
            continue
        if chrom == 'chrM':
            continue
        if line_dict['junction_letters'] == '-':
            PTES_logger.error('PE input, junction type -1 is present!')
            continue
        if abs(line_dict['donor_ss'] - line_dict['acceptor_ss']) > 1000000 \
                or chain == '+' and line_dict['donor_ss'] < line_dict['acceptor_ss'] \
                or chain == '-' and line_dict['donor_ss'] > line_dict['acceptor_ss']:
            if line_dict['junction_letters'] == 'GT/AG':
                mates_gtag['non-chim'] += 1
            else:
                mates_nc['non-chim'] += 1
            continue
        read_name = line_dict['read_name']
        mate1, mate2 = return_mates(cigar1=line_dict['cigar1'],
                                    coord1=line_dict['coord1'],
                                    cigar2=line_dict['cigar2'],
                                    coord2=line_dict['coord2'],
                                    chain=chain)
        interval_intersection = mate_intersection(mate1, mate2)
        if line_dict['junction_letters'] == 'GT/AG':
            mates_gtag[interval_intersection] += 1
        else:
            mates_nc[interval_intersection] += 1
        if interval_intersection == 'outside':
            outside_list.append(line)
        if line_dict['donor_ss'] in gtf_donors[chrom]:
            annot_donors += 1
        if line_dict['acceptor_ss'] in gtf_acceptors[chrom]:
            annot_acceptors += 1


PTES_logger.info('Reading STAR output... done')  
PTES_logger.info('Inside GT/AG: %i' % mates_gtag['inside'])
PTES_logger.info('Inside other: %i' % mates_nc['inside'])
PTES_logger.info('Outside GT/AG: %i' % mates_gtag['outside'])
PTES_logger.info('Outside other: %i' % mates_nc['outside'])
PTES_logger.info('Intron too large GT/AG: %i' % mates_gtag['non-chim'])
PTES_logger.info('Intron too large other: %i' % mates_nc['non-chim'])
PTES_logger.info('Annot donors: %i' % annot_donors)
PTES_logger.info('Annot acceptors: %i' % annot_acceptors)

writeln_to_file(''.join(outside_list), outside_name, folder = path_to_file)   