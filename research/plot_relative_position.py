# Takes bedtools intersect -wo output for features and containers
# Plots relative location of features in containers
# For feature [a, b] and container [x, y]:
# P = (a-x)/((y-x)-(b-a))
# p is float in [0,1], or -1 for non-intersecting intervals, or 2 for features larger than genes (just symmetry)


# Imports
import argparse
import os
import json

from interval import interval
import matplotlib.pyplot as plt
from collections import defaultdict
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
PTES_logger.info('Creating intersection file... ')
intersection_name = os.path.join(args.output, os.path.basename(args.features)+'.intersect')
cmd = 'bedtools intersect -a %s -b %s -s -wo > %s' % (args.features,
                                                      args.genes,
                                                      intersection_name,
                                                      )
shell_call(cmd)
PTES_logger.info('Creating intersection file... done')
PTES_logger.info('Reading intersection file... ')

p_dict = {}
gene_p_dict = defaultdict(list)
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
        gene_name = line_list[b_start+3]
        gene_tuple = (chrom1, line_list[b_start+5], line_list[b_start + 1], line_list[b_start + 2], gene_name)
        feature_len = get_interval_length(feature_interval)
        feature_len_list.append(feature_len)
        gene_len = get_interval_length(gene_interval)
        gene_len_list.append(gene_len)
        if feature_len <= gene_len:
            if args.strand:
                feature_list = [chrom1, line_list[5], line_list[1], line_list[2]]  # chrom, chain, start, end
                try:
                    container_strand = line_list[b_start + 5]
                    p = count_relative_position(
                        feature=feature_interval,
                        container=gene_interval,
                        container_strand=container_strand)
                    feature_tuple = tuple(feature_list + [p])
                    gene_p_dict[gene_tuple].append(feature_tuple)
                except IndexError:
                    PTES_logger.error('No strand found')
                    PTES_logger.error('BED6 format is required for plotting strand-specific position')
                    continue
            else:
                feature_list = [chrom1, line_list[1], line_list[2]]  # chrom, start, end
                p = count_relative_position(feature=feature_interval,
                                            container=gene_interval)
                feature_tuple = tuple(feature_list + [p])
                gene_p_dict[gene_tuple].append(feature_tuple)
        else:
            if args.strand:
                feature_tuple = (chrom1, line_list[5], line_list[1], line_list[2], 2)  # chrom, chain, start, end
            else:
                feature_tuple = (chrom1, line_list[1], line_list[2], 2)  # chrom, start, end

            gene_p_dict[gene_tuple].append(feature_tuple)


PTES_logger.info('Reading intersection file... done')
PTES_logger.info('Saving output to file... ')
PTES_logger.info('Length of dict is %i ' % len(gene_p_dict))
with open(os.path.join(args.output, args.prefix+'_gene_p_dict.json'), 'w') as gene_p_json:
    json.dump({str(k): v for k, v in gene_p_dict.items()}, gene_p_json, indent=2)
with open(os.path.join(args.output, args.prefix+'_lengths.txt'), 'w') as out_file:
    out_file.write('Feature lengths: \n')
    out_file.write(','.join(map(str, feature_len_list)) + '\n')
    out_file.write('Gene lengths: \n')
    out_file.write(','.join(map(str, gene_len_list)) + '\n')
PTES_logger.info('Saving output to file... done')

PTES_logger.info('Plotting... ')
all_values = []
for gene, values in gene_p_dict.items():
    for value in values:
        all_values.append(value[-1])

fig = plt.figure(figsize=(6,4))
plt.suptitle('Features: %s;' % args.features + '\n' + 'containers: %s' % args.genes)
ax1 = fig.add_subplot(111)
ax1.hist(all_values, bins=100)
ax1.set_xlim(0,1)
fig.subplots_adjust(top=0.9)


#fig.tight_layout()
plt.savefig(os.path.join(args.output, args.prefix+'_relative_position.png'))

PTES_logger.info('Plotting... done')




