# Takes paired-end STAR output Chimeric.out.junction
# Finds "mate-outside" and "mate-inside", with GT/AG and unique mapped reads,
# Counts x / (x + y) where x is N(mate-outside) and y is N(mate-inside)
# Copies lines with mates outside to mates_outside.junction to make bigBed files for UCSC genome browser

# Imports
from collections import defaultdict
import argparse
import random


import pandas as pd
import numpy as np
from interval import interval
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call, make_dir
from ptes.ptes import annot_junctions, \
    mate_intersection, get_read_interval, dict_to_interval, one_interval, \
    star_line_dict
from ptes.ucsc.ucsc import list_to_dict, get_track_list, make_bed_folder, to_bigbed


# Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str,
                    help="STAR output, Chimeric.out.junction \
                        OR list of such files")
parser.add_argument("-s", "--sam", type=str,
                    help="Filtered STAR SAM output, with read_names same as in Chimeric.out.junction OR list")
parser.add_argument("-o", "--output", type=str,
                    help="Output folder for results")
parser.add_argument("-l", "--list", type=str,
                    help="Enables list input mode. Options input, sam, tag MUST be lists")
parser.add_argument("-g", "--genome", type=str,
                    default='/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',
                    help="Absolute path to genome file")
parser.add_argument("-gtf", "--gtf_annot", type=str,
                    default='/home/sunnymouse/Human_ref/hg19_exons.gtf',
                    help="Absolute path to genome file")
parser.add_argument("-t", "--tag", type=str,
                    default='ENCODE',
                    help="Tag name for grouping results, i.e. ENCODE id OR list of tags")
args = parser.parse_args()

# Functions


def nonchim_input(sam_name, chim_name):
    """
    Reads STAR non-chimeric output
    To make it fast, skips lines with names that are absent in chimeric output 
    :param sam_name: STAR output Aligned.out.sam
    :param chim_name: STAR output Chimeric.out.junction.filtered
    :return: defaultdict sam_dict: key is read_name, value is list of dicts with SAM attrs 
    """
    shell_call('cat %s | cut -f 10 > %s.chim_names' % (chim_name, chim_name))
    with open('%s.chim_names' % chim_name, 'rb') as names_file:
        names_list = [x.strip(b'\n') for x in names_file.readlines()]
        names_set = set(names_list)  # only names with chimeric output

    sam_dict = defaultdict(list)
    with open(sam_name, 'rb') as input_file:
        for line in input_file:
            row = line.strip().split(b'\t')
            read_name = row[0]
            if read_name not in names_set:
                continue
            flag = int(row[1])
            if flag & 16 == 0:
                chain = '+'
            else:
                chain = '-'
            sam_attrs = {
                'chain': chain,
                'chrom': row[2],
                'leftpos': row[3],
                'cigar': row[5],
                'NH': int(row[11].lstrip(b'NH:i:')),
            }
            sam_dict[read_name].append(sam_attrs)
    return sam_dict


def chim_input(chim_name):
    """
    Reads STAR chimeric output
    :param chim_name: name of file Chimeric.out.junction.filtered
    :return: List of dicts ready for making a DataFrame,
    each dict is one read, mapped as mate-inside or mate-outside
    """
    annot_donors = 0
    annot_acceptors = 0
    mates = {'inside': 0, 'outside': 0, 'non-chim': 0}
    read_names_list = []

    with open(chim_name, 'rb') as input_file:
        for line in input_file:
            line_dict = star_line_dict(line=line)
            chrom = line_dict['chrom1']
            if chrom == 'chrM':
                continue
            chain = line_dict['chain1']
            read_name = line_dict['read_name']
            if line_dict['junction_letters'] == '-':
                PTES_logger.error('PE input, junction type -1 is present!')
                PTES_logger.warning('Skipping row...')
                continue
            if abs(line_dict['donor_ss'] - line_dict['acceptor_ss']) > 1000000 \
                    or chain == '+' and line_dict['donor_ss'] < line_dict['acceptor_ss'] \
                    or chain == '-' and line_dict['donor_ss'] > line_dict['acceptor_ss']:
                mates['non-chim'] += 1
                continue

            annot_donor = 0
            annot_acceptor = 0
            if line_dict['donor_ss'] in gtf_donors[chrom]:
                annot_donor = 1
                annot_donors += 1
            if line_dict['acceptor_ss'] in gtf_acceptors[chrom]:
                annot_acceptor = 1
                annot_acceptors += 1

            chim_part1 = get_read_interval(cigar=line_dict['cigar1'], leftpos=line_dict['coord1'])
            chim_part2 = get_read_interval(cigar=line_dict['cigar2'], leftpos=line_dict['coord2'])
            mate1 = one_interval(dict_to_interval(chim_part1) | dict_to_interval(chim_part2))
            for mate in sam_dict[read_name]:
                if mate['cigar'] == line_dict['cigar1'] or mate['cigar'] == line_dict['cigar2']:
                    continue
                if mate['NH'] > 1:
                    nh_chroms = 0  # check if mapping to this chromosome is unique
                    for mapping in sam_dict[read_name]:
                        if mapping['chrom'] == chrom:
                            nh_chroms += 1
                    if nh_chroms > 1:
                        continue
                if mate['chrom'] == chrom and mate['chain'] != chain:
                    mate2_dict = get_read_interval(cigar=mate['cigar'], leftpos=mate['leftpos'])
                    mate2 = one_interval(dict_to_interval(mate2_dict))
                    interval_intersection = mate_intersection(mate1, mate2)
                    mates[interval_intersection] += 1
                    junc_dict[(chrom,
                               chain,
                               line_dict['donor_ss'],
                               line_dict['acceptor_ss'],
                               )].append((chim_part1, chim_part2, mate2_dict, interval_intersection))
                    if interval_intersection == 'outside':
                        mate_dist = int(min(
                            abs(mate1[0].inf - mate2[0].sup),
                            abs(mate2[0].inf - mate1[0].sup),
                        ))
                    else:
                        mate_dist = 0
                    read_attrs = {
                        'read_name': read_name,
                        'chain': chain,  # chain of chimeric junction
                        'chrom': chrom,
                        'donor': line_dict['donor_ss'],
                        'acceptor': line_dict['acceptor_ss'],
                        'annot_donor': annot_donor,
                        'annot_acceptor': annot_acceptor,
                        'consensus': line_dict['junction_letters'],
                        'chim_dist': abs(line_dict['donor_ss'] - line_dict['acceptor_ss']),
                        'mate_dist': mate_dist,
                        'type': interval_intersection,
                    }
                    read_names_list.append(read_attrs)

    PTES_logger.info('Inside: %i' % mates['inside'])
    PTES_logger.info('Outside: %i' % mates['outside'])
    PTES_logger.info('Intron too large: %i' % mates['non-chim'])
    PTES_logger.info('Annot donors: %i' % annot_donors)
    PTES_logger.info('Annot acceptors: %i' % annot_acceptors)
    return read_names_list


def digit_code(number, threshold=6):
    """
    Makes code of length 6, adds 0 to the left
    Example: 125 -> 000125
    :param number: integer number w. 6 or less digits
    :param threshold: number of digits, default 6
    :return: string of same number with 6 digits
    """
    current_len = len(str(number))
    if current_len >= threshold:
        return str(number)
    else:
        n_zeros = threshold-current_len
        return '0'*n_zeros+str(number)

# Main


# Exons GTF to junctions dict
PTES_logger.info('Reading GTF...')
gtf_exons_name = args.gtf_annot
gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)
PTES_logger.info('Reading GTF... done')

# non-iterative
make_dir(args.output)
path_to_file = args.output.rstrip('/')
junc_dict = defaultdict(list)

if args.list:
    with open(args.input, 'r') as chim_names_file:
        chim_names_list = [x.strip('\n') for x in chim_names_file.readlines()]
    with open(args.sam, 'r') as sam_names_file:
        sam_names_list = [x.strip('\n') for x in sam_names_file.readlines()]
    with open(args.tag, 'r') as tag_names_file:
        tag_names_list = [x.strip('\n') for x in tag_names_file.readlines()]
    pairs = zip(chim_names_list, sam_names_list, tag_names_list)
    PTES_logger.info('Enabled list mode')
    chim_junc_df_list = []

    for chim_name, sam_name, tag in pairs:
        PTES_logger.info('Input file %s' % chim_name)
        # Reading filtered STAR non-chim output
        PTES_logger.info('Reading STAR non-chimeric output...')
        sam_dict = nonchim_input(sam_name=sam_name, chim_name=chim_name)
        PTES_logger.info('Reading STAR non-chimeric output... done')
        # Reading filtered STAR output
        PTES_logger.info('Reading STAR chimeric output...')
        read_names_list = chim_input(chim_name=chim_name)
        PTES_logger.info('Reading STAR chimeric output... done')
        chim_df = pd.DataFrame(read_names_list)
        try:
            chim_df = chim_df[['read_name', 'chrom', 'chain',
                                'donor', 'acceptor', 'annot_donor',
                                'annot_acceptor', 'consensus',
                                'chim_dist', 'mate_dist',
                                'type']].sort_values(by='read_name').reset_index(drop=True)
            chim_df['id'] = tag
            chim_junc_df_list.append(chim_df)
        except KeyError:
            PTES_logger.warning('Creating reads dataframe... empty dataframe')

    chim_junc_df = pd.concat(chim_junc_df_list).reset_index(drop=True)

if not args.list:
    PTES_logger.info('Input file %s' % args.input)

    # Reading filtered STAR non-chim output
    PTES_logger.info('Reading STAR non-chimeric output...')
    sam_dict = nonchim_input(sam_name=args.sam, chim_name=args.input)
    PTES_logger.info('Reading STAR non-chimeric output... done')

    # Reading filtered STAR output
    PTES_logger.info('Reading STAR chimeric output...')

    read_names_list = chim_input(chim_name=args.input)

    PTES_logger.info('Reading STAR chimeric output... done')

    chim_junc_df = pd.DataFrame(read_names_list)
    try:
        chim_junc_df = chim_junc_df[[
            'read_name', 'chrom', 'chain',
            'donor', 'acceptor', 'annot_donor',
            'annot_acceptor', 'consensus',
            'chim_dist', 'mate_dist', 'type',
        ]].sort_values(by='read_name').reset_index(drop=True)
        chim_junc_df['id'] = args.tag
    except KeyError:
        PTES_logger.warning('Creating reads dataframe... empty dataframe')


# Creating reads dataframe
PTES_logger.info('Creating reads dataframe...')
chim_junc_df.to_csv('%s/chim_junc_df.csv' % path_to_file, sep='\t')
filtering_rule = 'consensus=="GT/AG" & chim_dist < 10000 & mate_dist < 10000'
df_new = chim_junc_df.query(filtering_rule)
df_new['annot'] = df_new.annot_acceptor + df_new.annot_donor
z = df_new.groupby(['chrom', 'chain', 'donor', 'acceptor', 'type', 'annot']).size().reset_index(name='counts')
zz = z.pivot_table(index=['chrom', 'chain', 'donor', 'acceptor'], columns=['type'], values='counts', fill_value=0)
if 'outside' in zz:
    zz['all'] = zz.inside + zz.outside
else:
    zz['outside'] = 0
    zz['all'] = zz.inside
zz = zz.sort_values(by='all', ascending=False)
zz.to_csv('%s/chim_types.csv' % path_to_file, sep='\t')


PTES_logger.info('Creating reads dataframe... done')
PTES_logger.info('Making BED files...')
if args.list:
    bed_prefix = 'ENCODE'+str(len(pairs))
else:
    bed_prefix = args.tag
bed_name = '%s.bed' % bed_prefix  # only track lines
coord_name = '%s.coords.csv' % bed_prefix  # table with windows to paste into GB and with descriptions
info_name = '%s.track' % bed_prefix  # file to submit to GB

folder_name = '%s/bed/' % path_to_file
make_bed_folder(folder_name=folder_name,
                bed_name=bed_name,
                coord_name=coord_name,
                info_name=info_name,
                data_desc=bed_prefix,
                )

writeln_to_file('\t'.join(
    ['#window', 'inside', 'outside', 'annot', 'n_samples','codes']), coord_name, folder=folder_name)

single_name = '%s.single.bed' % bed_prefix   # one line - both mates, for intersecting
init_file(filename=single_name, folder=folder_name)

index_list = ['chrom', 'chain', 'donor', 'acceptor']

annot_table = z.pivot_table(index=index_list, values='annot')  # info about annotation
id_groups = df_new.groupby(index_list).apply(lambda x: x.id.nunique()).reset_index(name='id_counts')
id_table = id_groups.pivot_table(index=index_list, values='id_counts')  # info about ids
num = 0
for key, value in zz.iterrows():   # unique chimeric junctions
    chrom = key[0]  # key is (chrom, chain, donor_ss, acceptor_ss)
    chain = key[1]
    windows_min = []
    windows_max = []
    codes = []
    for chim1, chim2, nonchim, mate_type in junc_dict[key]:  # unique reads w. this junction
        num += 1
        code = digit_code(number=num)   # every unique number will be 6-digit
        codes.append(code)
        bed1 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=chim1,
                              name='%s_%s_chim1' % (code, mate_type),
                              color='r')
        bed2 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=chim2,
                              name='%s_%s_chim2' % (code, mate_type),
                              color='r')
        bed3 = get_track_list(chrom=chrom,
                              chain=chain,
                              read_dict=nonchim,
                              name='%s_%s_mate2' % (code, mate_type),
                              color='b')
        # Making BED file with one row for pair of mates
        interval_list = map(lambda x: dict_to_interval(x, put_n=False), [chim1, chim2, nonchim])
        single_interval = interval()
        for part in interval_list:
            single_interval = single_interval | part
        single_interval_list = [y for y in single_interval.components]
        single_track = get_track_list(
            chrom=chrom,
            chain=chain,
            read_dict=list_to_dict(single_interval_list),
            name='%s_%s' % (code, mate_type),
            color='255,0,255')   # for checking in GB that intervals are same
        writeln_to_file('\t'.join(single_track), single_name, folder=folder_name)

        for track_list in [bed1, bed2, bed3]:
            windows_min.append(int(track_list[1]))  # track_list[1] is chromStart, track_list[2] is chromEnd
            windows_max.append(int(track_list[2]))
            writeln_to_file('\t'.join(track_list), bed_name, folder=folder_name)

    window = (chrom,   # one window for junction
              min(windows_min) - 200,
              max(windows_max) + 200)
    description = '\t'.join(map(str,[   # one description for junction
        value.inside,
        value.outside,
        annot_table.loc[key].annot,
        id_table.loc[key].id_counts,
        '%s-%s' % (codes[0], codes[-1]),
                                    ]
                                )
                            )
    writeln_to_file('%s:%i-%i\t' % window + description, coord_name, folder=folder_name)

PTES_logger.info('Making BED files... done')
to_bigbed(bed_name=bed_name, folder_name=folder_name)
to_bigbed(bed_name=single_name, folder_name=folder_name)








