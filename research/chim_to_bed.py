# Takes STAR Chimeric.out.junction.filtered output
# Makes BED file for each chimeric junction: drawing ordinary introns

# Imports
from collections import defaultdict
import argparse
import random
import os
import errno

import pybedtools

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call
from ptes.ptes import get_read_interval, dict_to_interval, one_interval, \
    interval_to_string, get_interval_length
from ptes.ucsc.ucsc import list_to_dict, get_track_list, make_bed_folder, to_bigbed

### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="STAR output, Chimeric.out.junction")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-p", "--panhandles", type=str,
                    default='/home/sunnymouse/projects/PTES/panhandles/pairs.ucsc.bed',
                    help="BED file with panhandles")
parser.add_argument("-gtf", "--gtf_annot", type=str,
                    default='/home/sunnymouse/Human_ref/hg19_genes.gtf',
                    help="Absolute path to genome file")
parser.add_argument("-t", "--tag", type=str,
                    default='ENCODE',
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions

# Main

PTES_logger.info('Input file %s' % args.input)
# Reading filtered STAR output
PTES_logger.info('Reading STAR chimeric output...')
input_name = args.input
try:
    os.makedirs(args.output)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass
path_to_file = args.output.rstrip('/')

bed_name = '%s.bed' % args.tag  # only track lines
coord_name = '%s.coords.csv' % args.tag  # table with windows to paste into GB and with descriptions
info_name = '%s.track' % args.tag  # file to submit to GB
folder_name = '%s/bed/' % path_to_file
make_bed_folder(folder_name=folder_name,
                bed_name=bed_name,
                coord_name=coord_name,
                info_name=info_name,
                data_desc=args.tag)
writeln_to_file('\t'.join(['#window', 'donor', 'acceptor', 'chain', 'cigar1', 'cigar2']), coord_name, folder=folder_name)

with open(input_name, 'r') as input_file:
    for line in input_file:
        line_list = line.strip().split('\t')
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
        if chain == '+':
            if donor_ss < acceptor_ss or abs(donor_ss - acceptor_ss) > 1000000:
                continue
        elif chain == '-':
            if donor_ss > acceptor_ss or abs(donor_ss - acceptor_ss) > 1000000:
                continue
        if junction_type == '1':
            junction_letters = 'GT/AG'
        else:
            PTES_logger.warning('Not GT/AG')
            continue
        donor_part = dict_to_interval(get_read_interval(cigar=cigar1, leftpos=coord1), put_n=False)
        # interval of MD-intervals
        acceptor_part = dict_to_interval(get_read_interval(cigar=cigar2, leftpos=coord2), put_n=False)
        parts = donor_part | acceptor_part  # one interval of intervals
        parts_list = [x for x in parts.components]
        track_list = get_track_list(chrom=chrom,
                             chain=chain,
                             read_dict=list_to_dict(parts_list),
                             name='%i-%i' % (donor_ss, acceptor_ss),
                             )
        window = (chrom,
                  int(track_list[1]) - 200,
                  int(track_list[2]) + 200)
        writeln_to_file('\t'.join(track_list), bed_name, folder=folder_name)
        description = '\t'.join(map(str, [donor_ss, acceptor_ss, chain, cigar1, cigar2]))
        writeln_to_file('%s:%i-%i\t' % window + description, coord_name, folder=folder_name)

to_bigbed(bed_name=bed_name, folder_name=folder_name)
PTES_logger.info('Reading STAR chimeric output... done')
PTES_logger.info('Intersecting with panhandles...')
circles = pybedtools.BedTool('%s/%s' % (folder_name, bed_name))
panhandles = pybedtools.BedTool(args.panhandles).set_chromsizes('hg19')
genes = pybedtools.BedTool(args.gtf_annot)
#intersection = panhandles.intersect(circles)
#intersection_c = panhandles.intersect(circles, c=True)
results_dict = panhandles.randomstats(circles,
                                      iterations=10000,
                                      include_distribution=False,
                                      shuffle_kwargs={'chrom': True, 'genome': "hg19", 'incl': args.gtf_annot},
                                      debug=False,
                                      processes=10)

keys = ['self', 'other', 'actual', 'median randomized', 'normalized', 'percentile',
        'frac randomized above actual',
        'frac randomized below actual',
        'upper_97.5th',
        'lower_2.5th',]
with open('%s/results_dict' % folder_name, 'w') as res:
    for key in keys:
        res.write('%s: %s\n' % (key, results_dict[key]))
PTES_logger.info('Intersecting with panhandles... done')





