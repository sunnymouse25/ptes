# Takes bedtools intersect -wo output for features and containers
# Plots relative location of features in containers
# For feature [a, b] and container [x, y]:
# P = (a-x)/((y-x)-(b-a))
# p is float in [0,1], or -1 for non-intersecting intervals, or 2 for features larger than genes (just symmetry)


# Imports
import argparse
import os

from interval import interval
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import get_b_start, count_relative_position, get_interval_length


### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-p", "--prefix", type=str,
                    help="Output prefix for results")
parser.add_argument("-f", "--features", type=str,
                    help="Path to .BED file with features (small intervals)")
parser.add_argument("-g", "--genes", type=str,
                    help="Path to .BED file with genes (containers)")
parser.add_argument("-s", "--strand", type=str,
                    help="Enable strand-specific position mode, BED files should contain 6 or more fields")

args = parser.parse_args()

# Functions


# Main
make_dir(args.output)
path_to_file = args.output.rstrip('/')
PTES_logger.info('Creating intersection file... ')
intersection_name = '%s/%s' % (path_to_file, os.path.basename(args.features)+'.intersect')
cmd = 'bedtools intersect -a %s -b %s -s -wo > %s' % (args.features,
                                                   args.genes,
                                                   intersection_name,
                                                    )
shell_call(cmd)
PTES_logger.info('Creating intersection file... done')
PTES_logger.info('Reading intersection file... ')

p_dict = {}
feature_len_list = []
gene_len_list = []

with open(intersection_name, 'r') as intersect_file:
    for i, line in enumerate(intersect_file):
        line_list = line.strip().split()
        b_start = get_b_start(line)
        if not b_start:
            continue
        chrom1 = line_list[0]
        feature_interval = interval[int(line_list[1]), int(line_list[2])]
        gene_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]
        feature_len = get_interval_length(feature_interval)
        feature_len_list.append(feature_len)
        gene_len = get_interval_length(gene_interval)
        gene_len_list.append(gene_len)
        if feature_len <= gene_len:
            if args.strand:
                try:
                    container_strand = line_list[b_start + 5]
                    p_dict[feature_interval] = count_relative_position(feature=feature_interval,
                                                                       container=gene_interval,
                                                                       container_strand=container_strand)
                except IndexError:
                    PTES_logger.error('No strand found')
                    PTES_logger.error('BED6 format is required for plotting strand-specific position')
                    p_dict[feature_interval] = count_relative_position(feature=feature_interval,
                                                                       container=gene_interval,
                                                                       container_strand='.')
            else:
                p_dict[feature_interval] = count_relative_position(feature=feature_interval,
                                                                   container=gene_interval)
        else:
            p_dict[feature_interval] = 2


PTES_logger.info('Reading intersection file... done')

PTES_logger.info('Counting features that are absent in intersection file... ')

with open(args.features, 'r') as features_file:
    for i, line in enumerate(features_file):
        if line.startswith('#'):
            continue
        line_list = line.strip().split()
        chrom1 = line_list[0]
        feature_interval = interval[int(line_list[1]), int(line_list[2])]
        if not p_dict.get(feature_interval, None):
            p_dict[feature_interval] = -1

PTES_logger.info('Counting features that are absent in intersection file... done')

PTES_logger.info('Plotting... ')

fig = plt.figure(figsize=(8,12))
plt.suptitle('Features: %s;' % args.features + '\n' + 'containers: %s' % args.genes)
ax1 = fig.add_subplot(321)
ax1.hist(p_dict.values(), bins=150)
ax1.set(title='Relative location')

ax12 = fig.add_subplot(322)
ax12.boxplot(p_dict.values())
ax12.set(title='Relative location')

ax2 = fig.add_subplot(323)
ax2.hist(feature_len_list)
ax2.set_xlim(0, 10000)
ax2.set(title='Features length')

ax22 = fig.add_subplot(324)
ax22.boxplot(feature_len_list)
ax22.set_ylim(0, 200000)
ax22.set(title='Features length')

ax3 = fig.add_subplot(325)
ax3.hist(gene_len_list)
ax3.set_xlim(0, 500000)
ax3.set(title='Containers length')

ax32 = fig.add_subplot(326)
ax32.boxplot(gene_len_list)
ax3.set_ylim(0, 500000)
ax32.set(title='Containers length')

#fig.tight_layout()
plt.savefig('%s/%s_relative_position.png' % (path_to_file, args.prefix))

PTES_logger.info('Plotting... done')
PTES_logger.info('Saving output to file... ')
with open('%s/%s_relative_position.txt' % (path_to_file, args.prefix), 'w') as out_file:
    out_file.write('Relative positions: \n')
    out_file.write(','.join(map(str,p_dict.values())) + '\n')
    out_file.write('Feature lengths: \n')
    out_file.write(','.join(map(str, feature_len_list)) + '\n')
    out_file.write('Gene lengths: \n')
    out_file.write(','.join(map(str, gene_len_list)) + '\n')
PTES_logger.info('Saving output to file... done')



