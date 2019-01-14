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
import numpy as np
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import interval_to_string, randomize_interval, get_b_start

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


def return_category(a_interval, b_interval, logger=PTES_logger):
    """
    Makes A + B, returns category of A and B
    :param a_interval: type 'interval'
    :param b_interval: type 'interval'
    :param logger: logger
    :return: category
    """
    union = a_interval | b_interval
    intersection = a_interval & b_interval
    if intersection == interval():
        return 'no_overlap'
    elif union == a_interval:
        return 'b_in_a'
    elif union == b_interval:
        return 'a_in_b'
    elif a_interval in union and b_interval in union:
        return 'overlap'
    else:
        logger.warning('Unknown category')
        logger.warning('Input: A %s' % interval_to_string(a_interval))
        logger.warning('Input: B %s' % interval_to_string(b_interval))
        logger.warning('Intersection: %s' % interval_to_string(union))
        return 'unknown'


def parse_bedtools_wo(wo_outfile, logger=PTES_logger):
    """
    Bedtools -wo output to counts of intersections
    :param wo_outfile: bedtools intersect -wo output
    :param logger: logger
    :return: dict from pivot table: A_feature - N intersections for each category
    """
    df_list = []  # for lines to output
    with open(wo_outfile, 'r') as input_file:
        for line in input_file:
            line_list = line.strip().split()
            b_start = get_b_start(line, logger)
            a_interval = interval[int(line_list[1]), int(line_list[2])]
            b_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]
            category = return_category(a_interval=a_interval,
                                        b_interval=b_interval)
            output_list = line_list[:3] + [category]
            df_list.append(output_list)
    df = pd.DataFrame(df_list)
    df.columns = ['chrom', 'start', 'end', 'category']

    z = df.groupby(['chrom', 'start', 'end', 'category']).size().reset_index(name='counts')
    zz = z.pivot_table(index=['chrom', 'start', 'end'], columns=['category'], values='counts', fill_value=0)
    zz_dict = {k: v.tolist() for k, v in zz.iterrows()}
    return zz_dict


def random_shuffle_inside(intersect_name, n_iterations, random_folder, logger=PTES_logger):
    """
    Shuffles feature inside feature's gene
    :param intersect_name: bedtools intersect -wo file: feature with genes
    :param n_iterations: number of files to create
    :param random_folder: path to random subfolder
    :param logger: logger
    :return: files in random/ subfolder
    """
    n_list = np.empty((args.iterations, 0)).tolist()   # make list of 1000 empty lists

    with open(intersect_name, 'r') as intersect_file:
        for line in intersect_file:
            line_list = line.strip().split()
            b_start = get_b_start(line)
            if not b_start:
                continue
            chrom1 = line_list[0]
            b_interval = interval[int(line_list[1]), int(line_list[2])]
            gene_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]
            for n in range(n_iterations):
                random_interval = randomize_interval(small_i=b_interval, large_i=gene_interval)
                bed_list = [chrom1, int(random_interval[0].inf), int(random_interval[0].sup)]
                bed_line = '\t'.join(map(str,bed_list))
                n_list[n].append(bed_line)

        for n in range(n_iterations):
            out_name = random_folder + '/random_%i.bed' % n
            with open(out_name, 'w') as out_file:
                out_file.write('\n'.join(n_list[n]))


# Main


make_dir(args.output)
path_to_file = args.output.rstrip('/')
random_folder = path_to_file + '/random'
make_dir(random_folder)

category_name = 'categories.csv'
pivot_name = 'categories_pivot.csv'

PTES_logger.info('Reading B-feature file... ')
PTES_logger.info('B-feature file: %s ' % args.b_with_genes)

random_shuffle_inside(intersect_name=args.b_with_genes,
                      n_iterations=args.iterations,
                      random_folder=random_folder)

PTES_logger.info('Reading B-feature file... done')

PTES_logger.info('Reading input file... ')
PTES_logger.info('Input file: %s ' % args.input)
real_dict = parse_bedtools_wo(wo_outfile=args.input)

PTES_logger.info('Reading input file... done')

PTES_logger.info('Running random intersections... ')
random_dicts = []
for n in range(args.iterations):
    random_input = 'random_%i.bed' % n
    random_output = 'random_%i.bed.intersect' % n
    cmd = 'bedtools intersect -a %s -b %s/%s -wo > %s/%s' % (args.input,
                                                             random_folder,
                                                             random_input,
                                                             random_folder,
                                                             random_output)
    shell_call(cmd)
    random_dict = parse_bedtools_wo(wo_outfile=random_folder+'/'+random_output)
    random_dicts.append(random_dict)

PTES_logger.info('Running random intersections... done')

PTES_logger.info('Creating output files...')
# taking all keys together
key_set = set(real_dict.keys())
for random_dict in random_dicts:
    for k in random_dict.keys():
        key_set.add(k)

# concatenating all random dicts
all_random_data = defaultdict(list)

for key in key_set:
    for random_dict in random_dicts:
        if random_dict.get(key, None):
            all_random_data[key].append(random_dict[key])
        else:
            all_random_data[key].append([0,0,0])   # this feature has no overlaps with the current random file

# grouping results by categories, making output for histograms
all_random_groups = {}
mean_values_output = []

# summarizing number of intersections by keys
a_in_b_list = []
b_in_a_list = []
overlap_list = []

for i in range(args.iterations):
    a_in_b = 0
    b_in_a = 0
    overlap = 0
    for key in key_set:
        a_in_b += all_random_data[key][i][0]
        b_in_a += all_random_data[key][i][1]
        overlap += all_random_data[key][i][2]

    a_in_b_list.append(a_in_b)
    b_in_a_list.append(b_in_a)
    overlap_list.append(overlap)

real_a_in_b = 0
real_b_in_a = 0
real_overlap = 0
for key in real_dict:
    real_a_in_b += real_dict[key][0]
    real_b_in_a += real_dict[key][1]
    real_overlap += real_dict[key][2]

# plotting
fig = plt.figure(figsize=(8,12))
ax1 = fig.add_subplot(321)
ax1.hist(a_in_b_list)
ax1.axvline(real_a_in_b, color='r')
ax1.set(title='a_in_b');

ax12 = fig.add_subplot(322)
ax12.boxplot(a_in_b_list)
ax12.plot(1, real_a_in_b, color='w', marker='*', markeredgecolor='k')
ax12.set(title='a_in_b');

ax2 = fig.add_subplot(323)
ax2.hist(b_in_a_list)
ax2.axvline(real_b_in_a, color='r')
ax2.set(title='b_in_a');

ax22 = fig.add_subplot(324)
ax22.boxplot(b_in_a_list)
ax22.plot(1, real_b_in_a, color='w', marker='*', markeredgecolor='k')
ax22.set(title='b_in_a');

ax3 = fig.add_subplot(325)
ax3.hist(overlap_list)
ax3.axvline(real_overlap, color='r')
ax3.set(title='overlap');

ax32 = fig.add_subplot(326)
ax32.boxplot(overlap_list)
ax32.plot(1, real_overlap, color='w', marker='*', markeredgecolor='k')
ax32.set(title='overlap');

plt.savefig('%s/histograms.png' % path_to_file)

test_name = '%s/data_hist.test' % path_to_file
with open(test_name, 'w') as test_file:
    test_file.write('a_in_b' + '\n')
    test_file.write('real: %i' % real_a_in_b + '\n')
    test_file.write('random: '+','.join(map(str, a_in_b_list)) + '\n')

    test_file.write('b_in_a' + '\n')
    test_file.write('real: %i' % real_b_in_a + '\n')
    test_file.write('random: '+','.join(map(str, b_in_a_list)) + '\n')

    test_file.write('overlap' + '\n')
    test_file.write('real: %i' % real_overlap + '\n')
    test_file.write('random: '+','.join(map(str, overlap_list)) + '\n')

for key in key_set:
    all_random_groups[key] = zip(*all_random_data)
    # now we have: key - 3 lists 'a in b', 'b in a', 'overlap'
    # and real dict: key - 'a in b', 'b in a', 'overlap'
    # we need to sum by categories
    # also output for each key: mean number of intersections in every category
    means = list(np.around(np.mean(all_random_data[key], axis=0), decimals=3))
    out_key = '_'.join(map(str, list(key)))
    out_line = '\t'.join(map(str, [out_key] + real_dict.get(key, 0) + means))
    mean_values_output.append(out_line)

random_outname = 'categories_random.csv'
with open('%s/%s' % (path_to_file, random_outname), 'w') as random_outfile:
    random_outfile.write('\n'.join(mean_values_output))

PTES_logger.info('Creating output files... done')
PTES_logger.info('Remember to delete random subfolder')




