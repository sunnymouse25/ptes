# Takes bedtools intersect -wa -wb output,
# Attributes category for each row: 'a_in_b', 'b_in_a', 'overlap'
# Makes pivot 'unique_feature_A - n_a_in_b - n_b_in_a - n_overlap'
# Runs randomizing of B inside B's genes:
# bedtools intersect -a B_FEATURE -b GENES_NO_INTERSECTION -wo > output
# Plots histogram for number of random intersections with real value

# Imports
from collections import defaultdict
import argparse

from interval import interval
import pandas as pd
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import interval_to_string, get_b_start

### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="BED file, bedtools intersect -wa -wb output")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-iter", "--iterations", type=int,
                    default='1000',
                    help="Number of iterations for randomizing results")
parser.add_argument("-m", "--method", type=str,
                    nargs='+',
                    default=['inside', 'outside', 'bedtools', ],
                    help="Shuffling method(s): inside, outside, bedtools")

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
    col_list = zz.columns.tolist()
    for cat in ['a_in_b', 'b_in_a', 'overlap']:
        if cat not in col_list:
            zz[cat] = 0
    zz_dict = {k: v.tolist() for k, v in zz.iterrows()}
    return zz_dict

# Main


make_dir(args.output)
path_to_file = args.output.rstrip('/')
random_folder = path_to_file + '/random'

# category_name = 'categories.csv'
# pivot_name = 'categories_pivot.csv'

PTES_logger.info('Reading input file... ')
PTES_logger.info('Input file: %s ' % args.input)
real_dict = parse_bedtools_wo(wo_outfile=args.input)

PTES_logger.info('Reading input file... done')

PTES_logger.info('Running random intersections... ')
random_dicts = {}
for method in args.method:
    random_dicts[method] = []
    for n in range(args.iterations):
        random_input = '%s_%i.bed' % (method, n)
        random_output = '%s_%i.bed.intersect' % (method, n)
        cmd = 'bedtools intersect -a %s -b %s/%s -wo > %s/%s' % (args.input,
                                                                 random_folder,
                                                                 random_input,
                                                                 random_folder,
                                                                 random_output)
        shell_call(cmd)
        random_dict = parse_bedtools_wo(wo_outfile=random_folder+'/'+random_output)
        random_dicts[method].append(random_dict)

PTES_logger.info('Running random intersections... done')

PTES_logger.info('Creating output files...')

# taking all keys together
key_set = set(real_dict.keys())
for method in args.method:
    for random_dict in random_dicts[method]:
        for k in random_dict.keys():
            key_set.add(k)

# concatenating all random dicts
all_random_data = {}
for method in args.method:
    all_random_data[method] = defaultdict(list)
    for key in key_set:
        for random_dict in random_dicts[method]:
            if random_dict.get(key, None):
                all_random_data[method][key].append(random_dict[key])
            else:
                all_random_data[method][key].append([0,0,0]) # this feature has no overlaps with the current random file

# summarizing number of intersections by keys
a_in_b_list = {}
b_in_a_list = {}
overlap_list = {}

# grouping results by categories, making output for histograms
for method in args.method:
    a_in_b_list[method] = []
    b_in_a_list[method] = []
    overlap_list[method] = []
    for i in range(args.iterations):
        a_in_b = 0
        b_in_a = 0
        overlap = 0
        for key in key_set:
            a_in_b += all_random_data[method][key][i][0]
            b_in_a += all_random_data[method][key][i][1]
            overlap += all_random_data[method][key][i][2]

        a_in_b_list[method].append(a_in_b)
        b_in_a_list[method].append(b_in_a)
        overlap_list[method].append(overlap)

real_a_in_b = 0
real_b_in_a = 0
real_overlap = 0
for key in real_dict:
    real_a_in_b += real_dict[key][0]
    real_b_in_a += real_dict[key][1]
    real_overlap += real_dict[key][2]

# plotting
for method in args.method:
    fig = plt.figure(figsize=(8,12))
    plt.suptitle("Shuffling method %s" % method)
    ax1 = fig.add_subplot(321)
    ax1.hist(a_in_b_list[method])
    ax1.axvline(real_a_in_b, color='r')
    ax1.set(title='a_in_b');

    ax12 = fig.add_subplot(322)
    ax12.boxplot(a_in_b_list[method])
    ax12.plot(1, real_a_in_b, color='r', marker='*', markeredgecolor='k', markersize=15)
    ax12.set(title='a_in_b');

    ax2 = fig.add_subplot(323)
    ax2.hist(b_in_a_list[method])
    ax2.axvline(real_b_in_a, color='r')
    ax2.set(title='b_in_a');

    ax22 = fig.add_subplot(324)
    ax22.boxplot(b_in_a_list[method])
    ax22.plot(1, real_b_in_a, color='r', marker='*', markeredgecolor='k', markersize=15)
    ax22.set(title='b_in_a');

    ax3 = fig.add_subplot(325)
    ax3.hist(overlap_list[method])
    ax3.axvline(real_overlap, color='r')
    ax3.set(title='overlap');

    ax32 = fig.add_subplot(326)
    ax32.boxplot(overlap_list[method])
    ax32.plot(1, real_overlap, color='r', marker='*', markeredgecolor='k', markersize=15)
    ax32.set(title='overlap');

    plt.savefig('%s/histograms_%s.png' % (path_to_file, method))

# saving data for histograms
for method in args.method:
    hist_name = '%s/data_hist_%s' % (path_to_file, method)
    with open(hist_name, 'w') as hist_file:
        hist_file.write('a_in_b' + '\n')
        hist_file.write('real: %i' % real_a_in_b + '\n')
        hist_file.write('random: ' + ','.join(map(str, a_in_b_list[method])) + '\n')

        hist_file.write('b_in_a' + '\n')
        hist_file.write('real: %i' % real_b_in_a + '\n')
        hist_file.write('random: ' + ','.join(map(str, b_in_a_list[method])) + '\n')

        hist_file.write('overlap' + '\n')
        hist_file.write('real: %i' % real_overlap + '\n')
        hist_file.write('random: ' + ','.join(map(str, overlap_list[method])) + '\n')

'''
# for each key: mean number of intersections in every category
# really, no one needs this table
# if you want it, make sure to debug this code

all_random_groups = {}
mean_values_output = []
for method in args.method: 
    for key in key_set:
        all_random_groups[key] = zip(*all_random_data[method])
        # now we have: key - 3 lists 'a in b', 'b in a', 'overlap'
        # and real dict: key - 'a in b', 'b in a', 'overlap'
        # we need to sum by categories
        # also output for each key: mean number of intersections in every category
        means = list(np.around(np.mean(all_random_data[method][key], axis=0), decimals=3))
        out_key = '_'.join(map(str, list(key)))
        out_line = '\t'.join(map(str, [out_key] + real_dict.get(key, 0) + means))
        mean_values_output.append(out_line)
    
    random_outname = 'categories_random_%s.csv' % method
    with open('%s/%s' % (path_to_file, random_outname), 'w') as random_outfile:
        random_outfile.write('\n'.join(mean_values_output))
'''
PTES_logger.info('Creating output files... done')
PTES_logger.info('Remember to delete random subfolder')




