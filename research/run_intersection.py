# Takes single BED file
# Runs intersection with panhandles

# Imports
from collections import defaultdict
import argparse
import random

import pybedtools

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call


### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-a", "--a_input", type=str,
                    help="BED with single lines for each pair of mates")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-b", "--b_input", type=str,
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

PTES_logger.info('Input files: A %s ' % args.a_input)
PTES_logger.info('Input files: B %s ' % args.b_input)

circles = pybedtools.BedTool(args.a_input)
panhandles = pybedtools.BedTool(args.b_input).set_chromsizes('hg19')
genes = pybedtools.BedTool(args.gtf_annot)
#intersection = panhandles.intersect(circles)
#intersection_c = panhandles.intersect(circles, c=True)
PTES_logger.info('Intersecting circles with panhandles...')
ab_results_dict = panhandles.randomstats(circles,
                                      iterations=10000,
                                      include_distribution=False,
                                      shuffle_kwargs={'chrom': True, 'genome': "hg19", 'incl': args.gtf_annot},
                                      debug=False,
                                      processes=10)

PTES_logger.info('Intersecting circles with panhandles... done')
PTES_logger.info('Intersecting panhandles with circles... ')
ba_results_dict = circles.randomstats(panhandles,
                                      iterations=10000,
                                      include_distribution=False,
                                      shuffle_kwargs={'chrom': True, 'genome': "hg19", 'incl': args.gtf_annot},
                                      debug=False,
                                      processes=10)

PTES_logger.info('Intersecting panhandles with circles... done')
PTES_logger.info('Writing outputs...')

keys = ['self',
        'other',
        'actual',
        'median randomized',
        'normalized',
        'percentile',
        'frac randomized above actual',
        'frac randomized below actual',
        'upper_97.5th',
        'lower_2.5th',]
with open('results_dict', 'w') as res1, \
        open('results_dict', 'w') as res2:
    for key in keys:
        res1.write('%s: %s\n' % (key, ab_results_dict[key]))
        res2.write('%s: %s\n' % (key, ba_results_dict[key]))

PTES_logger.info('Writing outputs... done')


