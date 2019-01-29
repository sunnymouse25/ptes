# Takes bedtools intersect -wo output for features and containers
# Plots relative location of features in containers
# For feature [a, b] and container [x, y]:
# P = (a-x)/((y-x)-(b-a))
# p is float in [0,1], or -1 for non-intersecting intervals

# Imports
import argparse
import os

from interval import interval
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import get_b_start, count_relative_location


### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-f", "--features", type=str,
                    help="Path to .BED file with features (small intervals)")
parser.add_argument("-g", "--genes", type=str,
                    help="Path to .BED file with genes (containers)")

args = parser.parse_args()

# Functions


# Main
make_dir(args.output)
path_to_file = args.output.rstrip('/')
PTES_logger.info('Creating intersection file... ')
intersection_name = '%s/%s' % (path_to_file, os.path.basename(args.features)+'.intersect')
cmd = 'bedtools intersect -a %s -b %s -wo > %s' % (args.features,
                                                       args.genes,
                                                       intersection_name,
                                                       )
shell_call(cmd)
PTES_logger.info('Creating intersection file... done')
PTES_logger.info('Reading intersection file... ')

p_dict = {}
with open(intersection_name, 'r') as intersect_file:
    for i, line in enumerate(intersect_file):
        line_list = line.strip().split()
        b_start = get_b_start(line)
        if not b_start:
            continue
        chrom1 = line_list[0]
        feature_interval = interval[int(line_list[1]), int(line_list[2])]
        gene_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]
        p_dict[feature_interval] = count_relative_location(small_i=feature_interval,
                                                           large_i=gene_interval)


PTES_logger.info('Reading intersection file... done')

PTES_logger.info('Counting features that are absent in intersection file... ')

with open(args.features, 'r') as features_file:
    for i, line in enumerate(features_file):
        line_list = line.strip().split()
        chrom1 = line_list[0]
        feature_interval = interval[int(line_list[1]), int(line_list[2])]
        if not p_dict.get(feature_interval, None):
            p_dict[feature_interval] = -1

PTES_logger.info('Counting features that are absent in intersection file... done')

PTES_logger.info('Plotting... ')
fig = plt.figure(figsize=(4,6))
plt.suptitle('Features: %s;\n containers: %s' % (args.features, args.containers))
ax1 = fig.add_subplot(121)
ax1.hist(p_dict.values())

ax12 = fig.add_subplot(122)
ax12.boxplot(p_dict.values())

PTES_logger.info('Plotting... done')



