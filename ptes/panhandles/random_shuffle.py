# Suppose we have small interval (= feature) inside large interval (= container)
# Runs shuffling of intervals using 3 methods:
# 1. shuffle inside: randomly moves small interval inside large interval;
# 2. shuffle outside: randomly moves small interval inside the other large interval
# (approx. the same size OR coverage, same chromosome);
# 3. bedtools shuffle: randomly moves all features in the union of all containers
# TODO: test function choose_close

# Imports
from collections import defaultdict
import argparse
import random
import os

from interval import interval
import numpy as np

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import randomize_interval, get_b_start, get_interval_length, count_relative_position, one_interval

# Functions


def interval_to_bed_line(chrom, single_interval, name='.', score=0, strand='.'):
    """
    From interval to TAB-delimited line for .BED6 file
    :param chrom: chromosome
    :param single_interval: interval object with only one component
    :param name: name as in BED6
    :param score: int 0-1000
    :param strand: '+', '-' or '.'
    :return: TAB-delimited string with chrom, start, end, ame, score, strand
    """
    bed_interval = one_interval(I=single_interval) # just to be sure
    bed_list = [chrom, int(bed_interval[0].inf), int(bed_interval[0].sup), name, score, strand]
    bed_line = '\t'.join(map(str, bed_list))
    return bed_line


def choose_close(sorted_list, items='values', threshold=100):
    """
    Given sorted list of elements,
    for each element makes the list of all N close elements
    :param sorted_list: sorted list of elements
    :param items: list contains values or items as dict.items()
    :param threshold: 100 by default
    :return: dictionary: element - list of all close elements WITHOUT this element
    """
    interval_dict = {}  # gene - the other gene with approx. the same length
    len_list = len(sorted_list)
    for i, gene in enumerate(sorted_list):
        if i <= threshold:
            a = 0
            b = 2*threshold+1
        elif threshold < i < (len_list-threshold):
            a = i - threshold
            b = i + threshold + 1
        elif i >= (len_list-threshold):
            a = len_list - 2*threshold - 1
            b = len_list
        if items == 'values':
            list_slice = [x for x in sorted_list[a:b] if x != gene]
            interval_dict[gene] = list_slice
        if items == 'items':
            list_slice = [x[0] for x in sorted_list[a:b] if x != gene]
            interval_dict[gene[0]] = list_slice
    return interval_dict


# Main
def main():
    ### Arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--method", type=str,
                        nargs='+',
                        default=['inside', 'outside', 'bedtools', ],
                        help="Shuffling method(s): inside, outside, bedtools")
    parser.add_argument("-c", "--closest", type=str,
                        default='coverage',
                        help="Choose close elements for outside method by coverage or length")
    parser.add_argument("-o", "--output", type=str,
                        help="Output folder for results")
    parser.add_argument("-iter", "--iterations", type=int,
                        default='1000',
                        help="Number of iterations, default 1000")
    parser.add_argument("-f", "--features", type=str,
                        help="Path to .BED6 file with features (small intervals)")
    parser.add_argument("-g", "--genes", type=str,
                        help="Path to .BED6 file with genes (containers)")
    parser.add_argument("-s", "--chrom_sizes", type=str,
                        default='/home/sunnymouse/Human_ref/hg19.chrom.sizes',
                        help="The chrom_sizes file should be tab delimited and structured as follows: \
                        <chromName><TAB><chromSize>, use bedtools shuffle -h for details")

    args = parser.parse_args()
    make_dir(args.output)
    path_to_file = args.output.rstrip('/')
    random_folder = path_to_file + '/random'
    make_dir(random_folder)

    # Shuffling methods inside and outside:

    if 'inside' in args.method or 'outside' in args.method:
        if 'outside' in args.method:
            PTES_logger.info('Reading containers... ')
            PTES_logger.info('containers file: %s ' % args.genes)

            strand_dict = {}
            interval_dict = {}
            if args.closest == 'length':
                gene_dict = defaultdict(list)
                # open file with genes
                with open(args.genes, 'r') as genes_file:
                    for line in genes_file:
                        line_list = line.strip().split()
                        chrom = line_list[0]
                        gene_interval = interval[int(line_list[1]), int(line_list[2])]
                        gene_dict[chrom].append(gene_interval)
                        try:
                            strand = line_list[5]
                            strand_dict[gene_interval] = strand
                        except IndexError:
                            strand_dict[gene_interval] = '.'
                            PTES_logger.error('No strand found')
                            PTES_logger.error('BED6 format is required for choosing strand-specific position')

                # sort lists by gene length
                for key in gene_dict:  # key is chromosome
                    new_list = sorted(gene_dict[key], key=lambda x: get_interval_length(x))
                    interval_dict.update(choose_close(sorted_list=new_list))

            if args.closest == 'coverage':
                PTES_logger.info('Creating coverage file... ')

                cover_name = '%s/%s' % (random_folder, os.path.basename(args.features) + '.cov')
                cmd = 'bedtools coverage -a %s -b %s -s > %s' % (args.genes,
                                                                 args.features,
                                                                 cover_name,
                                                                 )
                shell_call(cmd)
                PTES_logger.info('Creating coverage file... done')

                cover_dict = {}
                with open(cover_name, 'r') as cover_file:
                    for line in cover_file:
                        line_list = line.strip().split()
                        gene_interval = interval[int(line_list[1]), int(line_list[2])]
                        cov = float(line_list[-1])
                        cover_dict[gene_interval] = cov
                        try:
                            strand = line_list[5]
                            strand_dict[gene_interval] = strand
                        except IndexError:
                            strand_dict[gene_interval] = '.'
                            PTES_logger.error('No strand found')
                            PTES_logger.error('BED6 format is required for choosing strand-specific position')

                new_list = sorted(cover_dict.items(), key=lambda x: x[1])
                interval_dict.update(choose_close(sorted_list=new_list, items='items'))

            PTES_logger.info('Reading containers... done')

        PTES_logger.info('Creating intersection file... ')

        intersection_name = '%s/%s' % (random_folder, os.path.basename(args.features)+'.intersect')
        cmd = 'bedtools intersect -a %s -b %s -wo -s > %s' % (args.features,
                                                              args.genes,
                                                              intersection_name,
                                                              )
        shell_call(cmd)
        PTES_logger.info('Creating intersection file... done')

        PTES_logger.info('Reading intersection file and shuffling... ')
        PTES_logger.info('intersection file: %s' % intersection_name)

        if 'inside' in args.method:
            n_list_inside = np.empty((args.iterations, 0)).tolist()  # make list of 1000 empty lists

        if 'outside' in args.method:
            n_list_outside = np.empty((args.iterations, 0)).tolist()  # make list of 1000 empty lists

        with open(intersection_name, 'r') as intersect_file:
            for i, line in enumerate(intersect_file):
                line_list = line.strip().split()
                b_start = get_b_start(line)
                if not b_start:
                    continue
                chrom1 = line_list[0]
                feature_interval = interval[int(line_list[1]), int(line_list[2])]
                gene_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]

                for n in range(args.iterations):
                    if 'inside' in args.method:
                        random_interval_inside = randomize_interval(small_i=feature_interval, large_i=gene_interval)
                        n_list_inside[n].append(interval_to_bed_line(chrom=chrom1,
                                                                     single_interval=random_interval_inside,
                                                                     name=line_list[3],
                                                                     strand=line_list[5]))

                    if 'outside' in args.method:
                        new_large_interval = random.choice(interval_dict[gene_interval])  # choose one of closest genes
                        new_strand = strand_dict[new_large_interval]
                        feature_len = get_interval_length(feature_interval)
                        gene_len = get_interval_length(gene_interval)
                        if feature_len <= gene_len:
                            try:
                                container_strand = line_list[b_start + 5]
                                relative_position = count_relative_position(feature=feature_interval,
                                                                            container=gene_interval,
                                                                            container_strand=container_strand)
                                random_interval_outside = randomize_interval(small_i=feature_interval,
                                                                             large_i=new_large_interval,
                                                                             large_i_strand=new_strand,
                                                                             same_position=True,
                                                                             p=relative_position)
                            except IndexError:
                                PTES_logger.error('No strand found')
                                PTES_logger.error('BED6 format is required for choosing strand-specific position')
                                relative_position = count_relative_position(feature=feature_interval,
                                                                            container=gene_interval)
                                random_interval_outside = randomize_interval(small_i=feature_interval,
                                                                             large_i= new_large_interval,
                                                                             same_position=True,
                                                                             p=relative_position)
                        else:
                            random_interval_outside = randomize_interval(small_i=feature_interval,
                                                                         large_i=new_large_interval,
                                                                         same_position=False,
                                                                         )
                        n_list_outside[n].append(interval_to_bed_line(chrom=chrom1,
                                                                      single_interval=random_interval_outside,
                                                                      name=line_list[3],
                                                                      strand=line_list[5]
                                                                      ))

        PTES_logger.info('Reading intersection file and shuffling... done')
        PTES_logger.info('Creating output files... ')
        for n in range(args.iterations):
            if 'inside' in args.method:
                out_name = random_folder + '/%s_%i.bed' % ('inside', n)
                with open(out_name, 'w') as out_file:
                    out_file.write('\n'.join(n_list_inside[n]))

            if 'outside' in args.method:
                out_name = random_folder + '/%s_%i.bed' % ('outside', n)
                with open(out_name, 'w') as out_file:
                    out_file.write('\n'.join(n_list_outside[n]))

    PTES_logger.info('Creating output files... done')
    # Shuffling method 3
    if 'bedtools' in args.method:
        PTES_logger.info('Running bedtools shuffle... ')
        for n in range(args.iterations):
            random_file = 'bedtools_%i.bed' % n
            cmd = 'bedtools shuffle -incl %s -i %s -g %s -chrom > %s/%s' % (args.genes,
                                                                            args.features,
                                                                            args.chrom_sizes,
                                                                            random_folder,
                                                                            random_file)
            shell_call(cmd)
        PTES_logger.info('Running bedtools shuffle... done')


if __name__ == "__main__":
    main()