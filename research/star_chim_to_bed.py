# Takes STAR Chimeric.out.junction output,
# converts into bed file,
# by default does not filter STAR output, there is an option to do it

# Imports
from collections import defaultdict
import argparse

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call, make_dir, digit_code
from ptes.ptes import annot_junctions, \
    mate_intersection, get_read_interval, dict_to_interval, one_interval, \
    star_line_dict
from ptes.ucsc.ucsc import list_to_dict, get_track_list, make_bed_folder, to_bigbed

### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="STAR Chimeric.out.junction output")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-f", "--filter", type=str,
                    help="Enable filtering for STAR output, creates .filtered file")
parser.add_argument("-bb", "--bigbed", type=str,
                    help="Create .bigBed file")
parser.add_argument("-t", "--tag", type=str,
                    default='ENCODE',
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions

# Main
make_dir(args.output)
path_to_file = args.output.rstrip('/')

if args.filter:
    filtered_name = args.input.strip() + '.filtered'
    shell_call("cat %s | awk '$1 ==$4 && $3 ==$6' > %s" % (args.input, filtered_name))
    input_name = filtered_name
else:
    input_name = args.input

PTES_logger.info('Input file: %s ' % input_name)

folder_name, bed_name, coord_name = make_bed_folder(
    prefix=args.tag,
    path_to_file=path_to_file)

with open(input_name, 'r') as input_file:
    for i, line in enumerate(input_file):
        line_dict = star_line_dict(line=line)
        if line_dict['chrom1'] == line_dict['chrom2'] \
                and line_dict['chain1'] == line_dict['chain2']:
            chrom = line_dict['chrom1']
            chain = line_dict['chain1']
        else:
            PTES_logger.error('Non-filtered STAR output')
            PTES_logger.error('Use awk "$1 ==$4 && $3 ==$6" to filter')
            continue
        if chrom == 'chrM':
            continue
        if line_dict['junction_letters'] == '-':
            PTES_logger.error('PE input, junction type -1 is present!')
            continue
        if abs(line_dict['donor_ss'] - line_dict['acceptor_ss']) > 1000000 \
                or chain == '+' and line_dict['donor_ss'] < line_dict['acceptor_ss'] \
                or chain == '-' and line_dict['donor_ss'] > line_dict['acceptor_ss']:
            continue
        read_name = line_dict['read_name']
        chim_part1 = get_read_interval(cigar=line_dict['cigar1'], leftpos=line_dict['coord1'])
        chim_part2 = get_read_interval(cigar=line_dict['cigar2'], leftpos=line_dict['coord2'])
        code = digit_code(number=i)  # every unique number will be 6-digit
        windows_min = []
        windows_max = []
        bed1 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=chim_part1,
                              name='%s_chim1' % code,
                              color='r')
        bed2 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=chim_part2,
                              name='%s_chim2' % code,
                              color='r')
        for track_list in [bed1, bed2]:
            windows_min.append(int(track_list[1]))  # track_list[1] is chromStart, track_list[2] is chromEnd
            windows_max.append(int(track_list[2]))
            writeln_to_file('\t'.join(track_list), bed_name, folder=folder_name)
        window = (chrom,  # one window for junction
                  min(windows_min) - 200,
                  max(windows_max) + 200)
        writeln_to_file('%s:%i-%i\t' % window + code, coord_name, folder=folder_name)

if args.bigbed:
    to_bigbed(bed_name=bed_name, folder_name=folder_name)
