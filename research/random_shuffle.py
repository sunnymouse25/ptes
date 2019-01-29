# Suppose we have small interval (= feature) inside large interval (= container)
# Runs shuffling of intervals using 3 methods:
# 1. shuffle inside: randomly moves small interval inside large interval;
# 2. shuffle outside: randomly moves small interval inside the other large interval
# (approx. the same size, same chromosome);
# 3. bedtools shuffle: randomly moves all features in the union of all containers

# Imports
from collections import defaultdict
import argparse
import random

from interval import interval
import numpy as np

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import interval_to_string, randomize_interval, get_b_start, get_interval_length

### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--method", type=str,
                    nargs='+',
                    default=['inside', 'outside', 'bedtools', ],
                    help="Shuffling method(s): inside, outside, bedtools")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-iter", "--iterations", type=int,
                    default='1000',
                    help="Number of iterations, default 1000")
parser.add_argument("-f", "--features", type=str,
                    help="Path to .BED file with features (small intervals)")
parser.add_argument("-g", "--genes", type=str,
                    help="Path to .BED file with genes (containers)")
parser.add_argument("-c", "--chrom_sizes", type=str,
                    default='/home/sunnymouse/Human_ref/hg19.chrom.sizes',
                    help="The chrom_sizes file should be tab delimited and structured as follows: \
                    <chromName><TAB><chromSize>, use bedtools shuffle -h for details")

args = parser.parse_args()

# Functions


def interval_to_bed_line(chrom, one_interval,):
    """
    From interval to TAB-delimited line for .BED file
    :param chrom: chromosome
    :param one_interval: interval object with only one component
    :return: TAB-delimited string with chrom, start, end
    """
    bed_list = [chrom, int(one_interval[0].inf), int(one_interval[0].sup)]
    bed_line = '\t'.join(map(str, bed_list))
    return bed_line


def choose_close(sorted_list, threshold=100):
    """
    Given sorted list of elements,
    for each element makes the list of all N close elements
    :param sorted_list: sorted list of elements
    :param threshold: 100 by default
    :return: dictionary: element - list of all close elements WITHOUT this element
    """
    interval_dict = {}  # gene - the other gene with approx. the same length
    len_list = len(sorted_list)
    for i, gene in enumerate(sorted_list):
        if i <= threshold:
            a = 0
            b = 2*threshold
        elif threshold < i < (len_list-threshold):
            a = i - threshold
            b = i + threshold + 1
        elif i >= (len_list-threshold):
            a = len_list - 2*threshold - 1
            b = len_list
        slice = [x for x in sorted_list[a:b] if x != gene]
        interval_dict[gene] = slice
    return interval_dict


# Main
make_dir(args.output)
path_to_file = args.output.rstrip('/')
random_folder = path_to_file + '/random'
make_dir(random_folder)

if 'outside' in args.method:
    PTES_logger.info('Reading containers... ')
    PTES_logger.info('containers file: %s ' % args.genes)

    gene_dict = defaultdict(list)
    # open file with genes
    with open(args.genes, 'r') as genes_file:
        for line in genes_file:
            line_list = line.strip().split()
            chrom = line_list[0]
            gene_interval = interval[int(line_list[1]), int(line_list[2])]
            gene_dict[chrom].append(gene_interval)

    # sort lists by gene length
    interval_dict = {}
    for key in gene_dict:   # key is chromosome
        new_list = sorted(gene_dict[key], key=lambda x: get_interval_length(x))
        interval_dict.update(choose_close(sorted_list=new_list))

    PTES_logger.info('Reading containers... done')

# Shuffling methods inside and outside:

if 'inside' in args.method or 'outside' in args.method:
    PTES_logger.info('Creating intersection file... ')

    intersection_name = args.features+'.intersect'
    cmd = 'bedtools intersect -a %s -b %s -wo > %s/%s' % (args.features,
                                                          args.genes,
                                                          path_to_file,
                                                          intersection_name,
                                                          )
    shell_call(cmd)
    PTES_logger.info('Creating intersection file... done')

    PTES_logger.info('Reading intersection file and shuffling... ')
    PTES_logger.info('intersection file: %s/%s' % (path_to_file, intersection_name))

    if 'inside' in args.method:
        n_list_inside = np.empty((args.iterations, 0)).tolist()  # make list of 1000 empty lists

    if 'outside' in args.method:
        n_list_outside = np.empty((args.iterations, 0)).tolist()  # make list of 1000 empty lists

    with open('%s/%s' % (path_to_file, intersection_name), 'r') as intersect_file:
        for i, line in enumerate(intersect_file):
            line_list = line.strip().split()
            b_start = get_b_start(line)
            if not b_start:
                continue
            chrom1 = line_list[0]
            b_interval = interval[int(line_list[1]), int(line_list[2])]
            gene_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]

            for n in range(args.iterations):
                if 'inside' in args.method:
                    random_interval_inside = randomize_interval(small_i=b_interval, large_i=gene_interval)
                    n_list_inside[n].append(interval_to_bed_line(chrom=chrom1, one_interval=random_interval_inside))

                if 'outside' in args.method:
                    random_interval_outside = randomize_interval(small_i=b_interval,
                                                                 large_i=random.choice(interval_dict[gene_interval]))
                    n_list_outside[n].append(interval_to_bed_line(chrom=chrom1, one_interval=random_interval_outside))
                    if i == 0 and n < 2:
                        print gene_interval
                        print random.choice(interval_dict[gene_interval])
                        print random_interval_outside
                        print interval_to_bed_line(chrom=chrom1, one_interval=random_interval_outside)

# outputting
    for n in range(args.iterations):
        if 'inside' in args.method:
            out_name = random_folder + '/%s_%i.bed' % ('inside', n)
            with open(out_name, 'w') as out_file:
                out_file.write('\n'.join(n_list_inside[n]))

        if 'outside' in args.method:
            out_name = random_folder + '/%s_%i.bed' % ('outside', n)
            with open(out_name, 'w') as out_file:
                out_file.write('\n'.join(n_list_outside[n]))

PTES_logger.info('Reading intersection file and shuffling... done')

# Shuffling method 3
PTES_logger.info('Running bedtools shuffle... ')

if 'bedtools' in args.method:
    for n in range(args.iterations):
        random_file = 'bedtools_%i.bed' % n
        cmd = 'bedtools shuffle -incl %s -i %s -g %s -chrom > %s/%s' % (args.genes,
                                                                        args.features,
                                                                        args.chrom_sizes,
                                                                        random_folder,
                                                                        random_file)
        shell_call(cmd)

PTES_logger.info('Running bedtools shuffle... done')

