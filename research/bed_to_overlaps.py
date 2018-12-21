# Takes bedtools intersect -wa -wb output,
# Attributes category for each row: 'a_in_b', 'b_in_a', 'overlap'
# Makes tables: 'feature A - category' and pivot 'unique_feature_A - n_a_in_b - n_b_in_a - n_overlap'
# Runs randomizing of B inside B's genes:
# bedtools intersect -a B_FEATURE -b GENES_NO_INTERSECTION -wa -wb output

# Imports
from collections import defaultdict
import argparse

from interval import interval
import pandas as pd

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call, make_dir
from ptes.ptes import interval_to_string, randomize_interval


### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="BED file, bedtools intersect -wa -wb output")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-bwg", "--b_with_genes", type=str,
                    help="bedtools intersect -a B_FEATURE -b GENES_NO_INTERSECTION -wa -wb output")
parser.add_argument("-iter", "--iterations", type=int,
                    default='1000',
                    help="Number of iterations for randomizing results")
parser.add_argument("-gtf", "--gtf_annot", type=str,
                    default='/home/sunnymouse/Human_ref/hg19_genes_coding_no_overlap.gtf',
                    help="Absolute path to genome file")
parser.add_argument("-t", "--tag", type=str,
                    default='ENCODE',
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions


def return_category(a_interval, b_interval):
    """
    Makes A + B, returns category of A and B
    :param a_interval: type 'interval'
    :param b_interval: type 'interval'
    :return: category
    """
    union = a_interval | b_interval
    if union == a_interval:
        return 'b_in_a'
    elif union == b_interval:
        return 'a_in_b'
    elif a_interval in union and b_interval in union:
        return 'overlap'
    else:
        PTES_logger.warning('Unknown category')
        PTES_logger.warning('Input: A %s' % interval_to_string(a_interval))
        PTES_logger.warning('Input: B %s' % interval_to_string(b_interval))
        PTES_logger.warning('Intersection: %s' % interval_to_string(union))
        return 'unknown'

# Main


PTES_logger.info('Input file: %s ' % args.input)

df_list = []  # for lines to output
b_dict = {}   # key is B-feature's (chrom,start,end), value is gene's interval


PTES_logger.info('Reading B-feature file... ')
with open(args.b_with_genes, 'r') as bwg_file:
    for line in bwg_file:
        line_list = line.strip().split('\t')
        chrom1 = line_list[0]
        start1 = line_list[1]
        end1 = line_list[2]
        start2 = int(line_list[4])
        end2 = int(line_list[5])
        b_dict[(chrom1, start1, end1)] = interval[start2, end2]
PTES_logger.info('Reading B-feature file... done')

PTES_logger.info('Reading input file... ')
with open(args.input, 'r') as input_file:
    for line in input_file:
        line_list = line.strip().split('\t')
        chrom1 = line_list[0]
        for i, elm in enumerate(line_list[3:], start=3):
            if elm == chrom1:
                b_start = i   # number of field where feature B starts
                break
        if not b_start:
            PTES_logger.error('Feature B start not found, skipping row...')
            continue
        a_interval = interval[int(line_list[1]), int(line_list[2])]
        b_interval = interval[int(line_list[b_start+1]), int(line_list[b_start+2])]
        real_category = return_category(a_interval=a_interval,
                                   b_interval=b_interval)
        b_key = (chrom1, line_list[b_start+1], line_list[b_start+2])
        gene_interval = b_dict[b_key]
        random_dict = {
            'a_in_b': 0,
            'b_in_a': 0,
            'overlap': 0,
            'unknown': 0,
        }
        for n in range(args.iterations):
            random_interval = randomize_interval(small_i=b_interval, large_i=gene_interval)
            random_category = return_category(a_interval=a_interval,
                                              b_interval=random_interval)
            random_dict[random_category] += 1
        output_list = line_list[:3] + \
                      [real_category, random_dict['a_in_b'], random_dict['b_in_a'], random_dict['overlap']]
        df_list.append(output_list)


PTES_logger.info('Reading input file... done')
PTES_logger.info('Creating output files...')

make_dir(args.output)
path_to_file = args.output.rstrip('/')
category_name = 'categories.csv'
pivot_name = 'categories_pivot.csv'

with open('%s/%s' % (path_to_file, category_name), 'w') as category_file:
    for track in df_list:
        category_file.write('\t'.join(map(str, track)) + '\n')

df = pd.DataFrame(df_list)
df.columns = ['chrom', 'start', 'end', 'real_category', 'rnd_a_in_b', 'rnd_b_in_a', 'rnd_overlap']

z = df.groupby(['chrom', 'start', 'end', 'real_category']).size().reset_index(name='counts')
zz = z.pivot_table(index=['chrom', 'start', 'end'], columns=['real_category'], values='counts', fill_value=0)
zz.to_csv('%s/%s' % (path_to_file, pivot_name), sep='\t')

PTES_logger.info('Creating output files... done')




