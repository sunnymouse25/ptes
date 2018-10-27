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
# from Bio import SeqIO
# from Bio.SeqRecord import SeqRecord

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call
from ptes.ptes import annot_junctions, \
    mate_intersection, get_read_interval, dict_to_interval, one_interval, \
    interval_to_string, get_interval_length
from ptes.ucsc.ucsc import list_to_dict, get_track_list, make_bed_folder, to_bigbed


### Arguments

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="STAR output, Chimeric.out.junction")
parser.add_argument("-s","--sam", type=str,
                    help="Filtered STAR SAM output, with read_names same as in Chimeric.out.junction")
parser.add_argument("-o","--output", type=str,
                    help="Output folder for results")
parser.add_argument("-g","--genome", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',
                    help="Absolute path to genome file")
parser.add_argument("-gtf","--gtf_annot", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/hg19_exons.gtf',
                    help="Absolute path to genome file")
parser.add_argument("-t","--tag", type=str,
                    default = 'ENCODE',
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions


# Exons GTF to junctions dict

PTES_logger.info('Reading GTF...')
gtf_exons_name = args.gtf_annot
gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)
PTES_logger.info('Reading GTF... done')

# Reading filtered STAR non-chim output
PTES_logger.info('Reading STAR non-chimeric output...')
shell_call('cat %s | cut -f 10 > %s.chim_names' % (args.input, args.input))
with open('%s.chim_names' % args.input, 'rb') as names_file:
    names_list = [x.strip(b'\n') for x in names_file.readlines()]
    names_set = set(names_list)

sam_name = args.sam
sam_dict = defaultdict(list)
with open(sam_name, 'rb') as input_file:
    for line in input_file:
        row = line.strip().split(b'\t')
        read_name = row[0]
        if read_name not in names_set:
            continue
        flag = int(row[1])
        chrom = row[2]
        leftpos = row[3]
        cigar = row[5]
        nh = int(row[11].lstrip(b'NH:i:'))
        if flag & 16 == 0:
            chain = '+'
        else:
            chain = '-'
        sam_attrs = {
            'chain': chain,
            'chrom' : chrom,
            'leftpos' : leftpos,
            'cigar' : cigar,
            'NH' : nh,
        }
        sam_dict[read_name].append(sam_attrs)
PTES_logger.info('Reading STAR non-chimeric output... done')

# Reading filtered STAR output
PTES_logger.info('Reading STAR chimeric output...')
input_name = args.input
path_to_file = args.output.rstrip('/')

annot_donors = 0
annot_acceptors = 0
mates = {'inside': 0, 'outside': 0, 'non-chim' : 0}
read_names_list = []
junc_dict = defaultdict(list)

with open(input_name, 'rb') as input_file:
    for line in input_file:
        line_list = line.strip().split(b'\t')
        chrom = line_list[0]
        donor_ss = int(line_list[1])    #donor splice site coord
        chain = line_list[2]
        acceptor_ss = int(line_list[4])    #acceptor splice site coord
        junction_type = line_list[6]   #junction type: -1=encompassing junction (between the mates), 1=GT/AG, 2=CT/AC
        read_name = line_list[9]
        coord1 = int(line_list[10])
        cigar1 = line_list[11]
        coord2 = int(line_list[12])
        cigar2 = line_list[13]
        if chain == '+':
            if donor_ss < acceptor_ss or abs(donor_ss - acceptor_ss) > 1000000:
                mates['non-chim'] += 1
                continue
        elif chain == '-':
            if donor_ss > acceptor_ss or abs(donor_ss - acceptor_ss) > 1000000:
                mates['non-chim'] += 1
                continue
        annot_donor = 0
        annot_acceptor = 0
        if donor_ss in gtf_donors[chrom]:
            annot_donor = 1
            annot_donors += 1
        if acceptor_ss in gtf_acceptors[chrom]:
            annot_acceptor = 1
            annot_acceptors += 1
        if junction_type == '1':
            junction_letters = 'GT/AG'
        elif junction_type == '2':
            junction_letters = 'CT/AC'
        else:
            junction_letters = '.'
        chim_part1 = get_read_interval(cigar=cigar1, leftpos=coord1)
        chim_part2 = get_read_interval(cigar=cigar2, leftpos=coord2)
        mate1 = one_interval(dict_to_interval(chim_part1) | dict_to_interval(chim_part2))
        mate2 = None
        for mate in sam_dict[read_name]:
            if mate['cigar'] == cigar1 or mate['cigar'] == cigar2:
                continue
            if mate['NH'] > 1:
                nh_chroms = 0   # check if mapping to this chromosome is unique
                for mapping in sam_dict[read_name]:
                    if mapping['chrom'] == chrom:
                        nh_chroms +=1
                if nh_chroms > 1:
                    continue
            if mate['chrom'] == chrom and mate['chain'] != chain:
                mate2_dict = get_read_interval(cigar=mate['cigar'], leftpos=mate['leftpos'])
                mate2 = one_interval(dict_to_interval(mate2_dict))
                interval_intersection = mate_intersection(mate1, mate2)
                mates[interval_intersection] += 1
                junc_dict[(chrom, chain, donor_ss, acceptor_ss)].append((chim_part1, chim_part2, mate2_dict))
                if interval_intersection == 'outside':
                    mate_dist = int(min(
                        abs(mate1[0].inf-mate2[0].sup),
                        abs(mate2[0].inf-mate1[0].sup),
                    ))
                else:
                    mate_dist = 0
                read_attrs = {
                    'read_name': read_name,
                    'chain': chain,   # chain of chimeric junction
                    'chrom' : chrom,
                    'donor' : donor_ss,
                    'acceptor' : acceptor_ss,
                    'annot_donor': annot_donor,
                    'annot_acceptor': annot_acceptor,
                    'consensus': junction_letters,
                    'chim_dist': abs(donor_ss-acceptor_ss),
                    'mate_dist': mate_dist,
                    'type' : interval_intersection,
                    }
                read_names_list.append(read_attrs)


PTES_logger.info('Reading STAR chimeric output... done')

print 'Inside: %i' % mates['inside']
print 'Outside: %i' % mates['outside']
print 'Intron too large: %i' % mates['non-chim']
print 'Annot donors: %i' % annot_donors
print 'Annot acceptors: %i' % annot_acceptors

PTES_logger.info('Creating reads dataframe...')

chim_junc_df = pd.DataFrame(read_names_list)
try:
    chim_junc_df = chim_junc_df[['read_name', 'chrom', 'chain',
                             'donor', 'acceptor', 'annot_donor',
                             'annot_acceptor', 'consensus',
                             'chim_dist', 'mate_dist', 'type']].sort_values(by='read_name').reset_index(drop=True)
    chim_junc_df.to_csv('%s/chim_junc_df.csv' % path_to_file, sep='\t')
    df_new = chim_junc_df.query('consensus=="GT/AG" & chim_dist < 10000 & mate_dist < 10000')
    z = df_new.groupby(['chrom', 'chain', 'donor', 'acceptor', 'type']).size().reset_index(name='counts')
    zz = z.pivot_table(index=['chrom', 'chain', 'donor', 'acceptor'], columns=['type'], values='counts', fill_value=0)
    zz['all'] = zz.inside + zz.outside
    zz = zz.sort_values(by='all', ascending=False)
    zz.to_csv('%s/chim_types.csv' % path_to_file, sep='\t')

    PTES_logger.info('Making BED files...')

    bed_name = '%s.bed' % args.tag  # only track lines
    coord_name = '%s.coords.csv' % args.tag  # table with windows to paste into GB and with descriptions
    info_name = '%s.track' % args.tag  # file to submit to GB
    folder_name = '%s/bed/' % path_to_file
    make_bed_folder(folder_name=folder_name,
                    bed_name=bed_name,
                    coord_name=coord_name,
                    info_name=info_name,
                    data_desc=args.tag)

    writeln_to_file('#window\tinside\toutside', coord_name, folder=folder_name)

    for key, value in zz.iterrows():
        chrom = key[0]
        chain = key[1]
        windows_min = []
        windows_max = []
        for chim1, chim2, nonchim in junc_dict[key]:
            num = random.randint(100,999)
            bed1 = get_track_list(chrom, chain, chim1, name='%i_chim_part1' % num, color='r')
            bed2 = get_track_list(chrom, chain, chim2, name='%i_chim_part2' % num, color='r')
            bed3 = get_track_list(chrom, chain, nonchim, name='%i_mate2' % num, color='b')
            for track_list in [bed1, bed2, bed3]:
                windows_min.append(int(track_list[1]))  # track_list[1] is chromStart, track_list[2] is chromEnd
                windows_max.append(int(track_list[2]))
                writeln_to_file('\t'.join(track_list), bed_name, folder=folder_name)
        window = (chrom,
                  min(windows_min) - 200,
                  max(windows_max) + 200)
        description = "%i\t%i" % (value.inside, value.outside)
        writeln_to_file('%s:%i-%i\t' % window + description, coord_name, folder=folder_name)

    to_bigbed(bed_name=bed_name, folder_name = folder_name)

    PTES_logger.info('Creating reads dataframe... done')
except KeyError:
    PTES_logger.warning('Creating reads dataframe... empty dataframe')




