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
from ptes.ucsc.ucsc import list_to_dict, get_track_list, make_bed_folder, to_bigbed, get_single_track

### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="STAR Chimeric.out.junction output")
parser.add_argument("-o", "--output", type=str,
                    help="Path for ./bed subfolder with results")
parser.add_argument("-f", "--filter", type=str,
                    help="Enable filtering for STAR output, creates .filtered file")
parser.add_argument("-s", "--sort", type=str,
                    help="Sort BED files")
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
    PTES_logger.info('Filtering STAR output...')
    filtered_name = args.input.strip() + '.filtered'
    shell_call("cat %s | awk '$1 ==$4 && $3 ==$6' > %s" % (args.input, filtered_name))
    input_name = filtered_name
    PTES_logger.info('Filtering STAR output... done')
else:
    input_name = args.input

PTES_logger.info('Input file: %s ' % input_name)

folder_name, bed_name, coord_name = make_bed_folder(
    prefix=args.tag,
    path_to_file=path_to_file)

junc_bed_name = '%s.junctions.bed' % args.tag  # for unique junctions
single_bed_name = '%s.single.bed' % args.tag   # single line for one chimeric junction

skipped = {'non-filtered': 0,
           'chrM': 0,
           'j_type-': 0,
           'non-chim': 0}
bed_list = []  # for outputting BED lines
coord_list = []   # for outputting coord lines
single_list = [] # for outputting BED lines, one row per one chimeric junction
junc_dict = {}

PTES_logger.info('Reading STAR output...')
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
            skipped['non-filtered'] +=1
            continue
        if chrom == 'chrM':
            skipped['chrM'] += 1
            continue
        if line_dict['junction_letters'] == '-':
            PTES_logger.error('PE input, junction type -1 is present!')
            skipped['j_type-'] += 1
            continue
        if abs(line_dict['donor_ss'] - line_dict['acceptor_ss']) > 1000000 \
                or chain == '+' and line_dict['donor_ss'] < line_dict['acceptor_ss'] \
                or chain == '-' and line_dict['donor_ss'] > line_dict['acceptor_ss']:
            skipped['non-chim'] += 1
            continue
        read_name = line_dict['read_name']
        chim_part1 = get_read_interval(cigar=line_dict['cigar1'], leftpos=line_dict['coord1'])
        chim_part2 = get_read_interval(cigar=line_dict['cigar2'], leftpos=line_dict['coord2'])
        junction_name = '%i_%i' % (line_dict['donor_ss'], line_dict['acceptor_ss'])   # name is chim junc coords
        code = digit_code(number=i)  # every unique number will be 6-digit
        windows_min = []
        windows_max = []
        bed1 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=chim_part1,
                              name='%s_%s_%s_chim1' % (chrom, junction_name, code),
                              color='r')
        bed2 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=chim_part2,
                              name='%s_%s_%s_chim2' % (chrom, junction_name, code),
                              color='r')
        if not junc_dict.get(junction_name, None):   # for each unique junction write the 1st representative chim pair
            junc_dict[junction_name] = []
            add_junc = True
        else:
            add_junc = False
        for track_list in [bed1, bed2]:
            windows_min.append(int(track_list[1]))  # track_list[1] is chromStart, track_list[2] is chromEnd
            windows_max.append(int(track_list[2]))
            bed_line = '\t'.join(track_list)+'\n'
            bed_list.append(bed_line)
            if add_junc:
                junc_dict[junction_name].append(bed_line)

        window = (chrom,
                  min(windows_min) - 200,
                  max(windows_max) + 200)
        coord_list.append('%s:%i-%i\t' % window + '\t'.join([junction_name, code]) + '\n')
        # Making BED file with one row for pair of mates
        single_track = get_single_track(read_dict_list=[chim_part1, chim_part2],
                                        kwargs={
                                        'chrom': chrom,
                                        'chain': chain,
                                        'name': '%s_%s_%s' % (chrom, junction_name, code),
                                        'color': '255,0,255'}  # for checking in GB that intervals are same
                                        )
        single_list.append('\t'.join(single_track)+'\n')

PTES_logger.info('Reading STAR output... done')
PTES_logger.info('Processed: %i rows' % i)
for key in skipped:
    PTES_logger.info('Skipped %s: %i rows' % (key, skipped[key]))
PTES_logger.info('Converted successfully: %i rows' % len(coord_list))

PTES_logger.info('Writing BED files...')
with open('%s/%s' % (folder_name, bed_name), 'w') as bed_file, \
        open('%s/%s' % (folder_name, coord_name), 'w') as coord_file, \
        open('%s/%s' % (folder_name, junc_bed_name), 'w') as junc_bed_file,\
        open('%s/%s' % (folder_name, single_bed_name), 'w') as single_bed_file:
    for track_list in bed_list:
        bed_file.write(track_list)
    for coord in coord_list:
        coord_file.write(coord)
    for key, value in junc_dict.items():
        junc_bed_file.write(value[0])
        junc_bed_file.write(value[1])
    for track_list in single_list:
        single_bed_file.write(track_list)

PTES_logger.info('Writing BED files... done')

if args.sort:
    PTES_logger.info('Sorting BED files...')
    for filename in [bed_name, junc_bed_name, single_bed_name]:
        shell_call('cat %s/%s | sort -k1,1 -k2,2n  > %s/%s.sorted' % (folder_name,
                                                                      filename,
                                                                      folder_name,
                                                                      filename))

    PTES_logger.info('Sorting BED files... done')

if args.bigbed:
    PTES_logger.info('Making bigBed...')
    to_bigbed(bed_name=bed_name, folder_name=folder_name)
    PTES_logger.info('Making bigBed... done')