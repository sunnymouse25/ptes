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
import matplotlib.pyplot as plt
import seaborn as sns
plt.switch_backend('agg')

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir
from ptes.ptes import interval_to_string, get_b_start

# Functions


def return_category(a_interval, b_interval, logger=PTES_logger):
    """
    Makes A + B, returns category of A and B
    :param a_interval: type 'interval'
    :param b_interval: type 'interval'
    :param logger: logger
    :return: category from categories = ['a_in_b', 'b_in_a', 'overlap', 'unknown', 'no_overlap']
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
    categories = ['a_in_b', 'b_in_a', 'overlap', 'unknown', 'no_overlap']
    df_dict = dict.fromkeys(categories, 0)  # for output
    with open(wo_outfile, 'r') as input_file:
        for line in input_file:
            line_list = line.strip().split()
            b_start = get_b_start(line, logger)
            if not b_start:
                continue
            a_interval = interval[int(line_list[1]), int(line_list[2])]
            b_interval = interval[int(line_list[b_start + 1]), int(line_list[b_start + 2])]
            category = return_category(a_interval=a_interval,
                                        b_interval=b_interval)
            df_dict[category] += 1
    return df_dict

# Main

def main():
    ### Arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="BED file, bedtools intersect -a A_FEATURE -b B_BEATURE -s -wo output")
    parser.add_argument("-f", "--features", type=str,
                        help="Path to .BED6 file with features (small intervals)")
    parser.add_argument("-o", "--output", type=str,
                        help="Output folder with random subfolder for results")
    parser.add_argument("-a", "--afeature", type=str,
                        help="Name of A_FEATURE, i.e. circles")
    parser.add_argument("-b", "--bfeature", type=str,
                        help="Name of B_BEATURE, i.e. panhandles")
    parser.add_argument("-iter", "--iterations", type=int,
                        default='1000',
                        help="Number of iterations for randomizing results")
    parser.add_argument("-m", "--method", type=str,
                        nargs='+',
                        default=['inside', 'outside', 'bedtools', ],
                        help="Shuffling method(s): inside, outside, bedtools")

    args = parser.parse_args()

    path_to_file = args.output.rstrip('/')
    random_folder = path_to_file + '/random'

    # category_name = 'categories.csv'
    # pivot_name = 'categories_pivot.csv'

    PTES_logger.info('Reading input file... ')
    PTES_logger.info('Input file: %s ' % args.input)
    real_dict = parse_bedtools_wo(wo_outfile=args.input)  # category - number of intersections

    PTES_logger.info('Reading input file... done')

    PTES_logger.info('Running random intersections... ')
    random_dicts = {}   # method - category - list of values
    categories = ['a_in_b', 'b_in_a', 'overlap']
    for method in args.method:
        random_dicts[method] = dict.fromkeys(categories)
        for k, _ in random_dicts[method].items():
            random_dicts[method][k] = []   # now we have separate empty lists as values
        for n in range(args.iterations):
            random_input = '%s/%s_%i.bed' % (random_folder, method, n)
            random_output = '%s/%s_%i.bed.intersect' % (random_folder, method, n)
            cmd = 'bedtools intersect -a %s -b %s -s -wo > %s' % (args.features,
                                                                  random_input,
                                                                  random_output)
            shell_call(cmd)
            random_dict = parse_bedtools_wo(wo_outfile=random_output)
            for cat in categories:
                random_dicts[method][cat].append(random_dict[cat])

    PTES_logger.info('Running random intersections... done')

    PTES_logger.info('Creating output files...')
    # plotting
    for method in args.method:
        fig = plt.figure(figsize=(12,6))
        plt.suptitle("Shuffling method %s" % method)
        ax1 = fig.add_subplot(131)
        sns.distplot(random_dicts[method]['a_in_b'], kde=False)
        ax1.axvline(real_dict['a_in_b'], color='r')
        ax1.set(title='%s_in_%s' % (args.afeature, args.bfeature));

        ax2 = fig.add_subplot(132)
        sns.distplot(random_dicts[method]['b_in_a'], kde=False)
        ax2.axvline(real_dict['b_in_a'], color='r')
        ax2.set(title='%s_in_%s' % (args.bfeature, args.afeature));

        ax3 = fig.add_subplot(133)
        sns.distplot(random_dicts[method]['overlap'], kde=False)
        ax3.axvline(real_dict['overlap'], color='r')
        ax3.set(title='overlap');

        plt.savefig('%s/histograms_%s.png' % (path_to_file, method))

    # saving data for histograms
    for method in args.method:
        hist_name = '%s/data_hist_%s' % (path_to_file, method)
        with open(hist_name, 'w') as hist_file:
            hist_file.write('%s_in_%s' % (args.afeature, args.bfeature) + '\n')
            hist_file.write('real: %i' % real_dict['a_in_b'] + '\n')
            hist_file.write('random: ' + ','.join(map(str, random_dicts[method]['a_in_b'])) + '\n')

            hist_file.write('%s_in_%s' % (args.bfeature, args.afeature) + '\n')
            hist_file.write('real: %i' % real_dict['b_in_a'] + '\n')
            hist_file.write('random: ' + ','.join(map(str, random_dicts[method]['b_in_a'])) + '\n')

            hist_file.write('overlap' + '\n')
            hist_file.write('real: %i' % real_dict['overlap'] + '\n')
            hist_file.write('random: ' + ','.join(map(str, random_dicts[method]['overlap'])) + '\n')

    PTES_logger.info('Creating output files... done')
    PTES_logger.info('Remember to delete random subfolder')


if __name__ == "__main__":
    main()

