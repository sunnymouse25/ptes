# Takes single-end STAR output mate1_Chimeric.out.junction
# Finds "mate-outside" and "mate-inside", with GT/AG and unique mapped reads,
# Creates reads table, junctions table, junc_dict with read intervals as json file

# Imports
from collections import defaultdict
import argparse
import os
import json
import gzip

import pandas as pd

from ptes.constants import PTES_logger
from ptes.lib.general import make_dir
from ptes.ptes import annot_junctions, \
    mate_intersection, get_read_interval, dict_to_interval, one_interval, \
    star_line_dict

# Functions


def sam_input(sam_name, chim_name, norm_dict):
    """
    Reads STAR non-chimeric output
    To make it faster, skips lines with names that are absent in chimeric output
    :param sam_name: STAR output Aligned.out.sam
    :param chim_name: STAR output Chimeric.out.junction.filtered
    :param norm_dict: defaultdict: key is (chrom,chain,donor,acceptor), value is dict: read_name: read_intervals
    :return: defaultdict sam_dict: key is read_name, value is list of dicts with SAM attrs
    """
    with open(chim_name, 'r') as chim_file:
        names_list = [x.strip('\n').split('\t')[9] for x in chim_file.readlines()]
        names_set = set(names_list)  # only names with chimeric output

    sam_dict = defaultdict(list)   # for searching for non-chimeric mates
    with open(sam_name, 'r') as input_file:
        for line in input_file:
            row = line.strip().split('\t')
            read_name = row[0]
            if read_name in names_set:
                sam_attrs = parse_sam_row(row)
                if sam_attrs:
                    sam_dict[read_name].append(sam_attrs)
                    if 'N' in sam_attrs['cigar']:   # read mapped with intron
                        read_dict = get_read_interval(cigar=sam_attrs['cigar'],
                                                       leftpos = sam_attrs['leftpos'],
                                                       output='dict')
                        if sam_attrs['chain'] == '+':
                            donor_ss = int(read_dict['N1'][0].inf - 1)   # counts first N as intron
                            acceptor_ss = int(read_dict['N1'][0].sup + 1)
                        elif sam_attrs['chain'] == '-':
                            donor_ss = int(read_dict['N1'][0].sup + 1)
                            acceptor_ss = int(read_dict['N1'][0].inf - 1)
                        norm_dict[(sam_attrs['chrom'],
                                   sam_attrs['chain'],
                                   donor_ss,
                                   acceptor_ss)][read_name].append(read_dict)

    return sam_dict


def parse_sam_row(row):
    """
    Takes line of SAM file and makes a dict of attributes
    :param row: list (single row splitted by TAB) of STAR output Aligned.out.sam
    :return: dict if not chrM else None
    """
    flag = int(row[1])
    if flag & 16 == 0:
        chain = '+'
    else:
        chain = '-'
    if row[2] == 'chrM':
        return None
    sam_attrs = {
        'chain': chain,
        'chrom': row[2],
        'leftpos': row[3],
        'cigar': row[5],
        'NH': int(row[11].lstrip('NH:i:')),
    }
    return sam_attrs


def chim_input(chim_name, gtf_donors, gtf_acceptors, sam_dict, junc_dict):
    """
    Reads STAR chimeric output
    :param chim_name: name of file Chimeric.out.junction.filtered
    :param gtf_donors: dict: chrom - set of donors coordinates (integers)
    :param gtf_acceptors: dict: chrom - set of acceptors coordinates (integers)
    :param sam_dict: output of sam_input function
    :param junc_dict: defaultdict, global, with interval dicts for BED files
    :return: List of dicts ready for making a DataFrame,
    each dict is one read, mapped as mate-inside or mate-outside
    """
    annot_donors = 0
    annot_acceptors = 0
    read_names_list = []
    skipped = {'non-filtered': 0,   # different chromosomes and/or chains
               'chrM': 0,      # mapping to chrM
               'j_type-': 0,   # junction between the mates, -1 in STAR output
               'non-chim': 0}   # STAR counts very long (>1Mb) junctions as chimeric

    with open(chim_name, 'r') as input_file:
        for i, line in enumerate(input_file):
            line_dict = star_line_dict(line=line)
            if line_dict['chrom1'] == line_dict['chrom2'] \
                    or line_dict['chain1'] == line_dict['chain2']:
                chrom = line_dict['chrom1']
                chain = line_dict['chain1']
            else:
                PTES_logger.error('Non-filtered STAR output')
                PTES_logger.error('Use awk "$1 ==$4 && $3 ==$6" to filter')
                skipped['non-filtered'] += 1
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
            annot_donor = 0
            annot_acceptor = 0
            if line_dict['donor_ss'] in gtf_donors[chrom]:
                annot_donor = 1
                annot_donors += 1
            if line_dict['acceptor_ss'] in gtf_acceptors[chrom]:
                annot_acceptor = 1
                annot_acceptors += 1

            mate_tuple = return_mate_tuple(line_dict=line_dict,
                                           second_mates=sam_dict[read_name],
                                           chrom=chrom,
                                           chain=chain)
            if mate_tuple:
                junc_dict[(chrom,
                           chain,
                           line_dict['donor_ss'],
                           line_dict['acceptor_ss'],
                           )].append({read_name: mate_tuple[:3]})   # read_intervals only
                interval_intersection = mate_tuple[3]
                if interval_intersection == 'outside':
                    mate_dist = count_mate_outside_dist(mate_tuple=mate_tuple,
                                                        chain=chain)
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
                    'letters_ss': line_dict['junction_letters'],
                    'chim_dist': abs(line_dict['donor_ss'] - line_dict['acceptor_ss']),
                    'mate_dist': mate_dist,
                    'type': interval_intersection,
                }
                read_names_list.append(read_attrs)
    PTES_logger.info('Processed: %i rows' % i)
    for key in skipped:
        PTES_logger.info('Skipped %s: %i rows' % (key, skipped[key]))
    PTES_logger.info('Converted successfully: %i rows' % len(read_names_list))
    PTES_logger.info('Annot donors: %i' % annot_donors)
    PTES_logger.info('Annot acceptors: %i' % annot_acceptors)

    return read_names_list


def return_mate_tuple(line_dict, second_mates, chrom, chain):
    """
    Takes single line from Chimeric.out.junction file
    :param line_dict: dictionary of elements in single line, output of star_line_dict
    :param second_mates: list of reads with the same name from SAM file
    :param chrom: chromosome
    :param chain: chain
    :return: None or tuple of 3 read intervals and type of mapping (mate_inside/mate_outside)
    """
    chim_part1 = get_read_interval(cigar=line_dict['cigar1'],
                                   leftpos=line_dict['coord1'],
                                   )
    chim_part2 = get_read_interval(cigar=line_dict['cigar2'],
                                   leftpos=line_dict['coord2'],
                                   )
    mate1 = one_interval(dict_to_interval(chim_part1) | dict_to_interval(chim_part2))
    for mate in second_mates:
        if mate['cigar'] == line_dict['cigar1'] \
                or mate['cigar'] == line_dict['cigar2']:
            continue
        if mate['NH'] > 1:
            nh_chroms = 0  # check if mapping to this chromosome is unique
            for mapping in second_mates:
                if mapping['chrom'] == chrom:
                    nh_chroms += 1
            if nh_chroms > 1:
                continue
        if mate['chrom'] == chrom and mate['chain'] != chain:
            mate2_dict = get_read_interval(cigar=mate['cigar'], leftpos=mate['leftpos'])
            mate2 = one_interval(dict_to_interval(mate2_dict))
            interval_intersection = mate_intersection(mate1, mate2)
            return chim_part1, chim_part2, mate2_dict, interval_intersection
    return None


def count_mate_outside_dist(mate_tuple, chain):
    """
    Counts distance between mates in "mate_outside", as minimal distance
    Configuration: "right" if non-chim mate is closer to 3'-end than chimeric mate,
    "left" if chimeric mate is closer to 3'-end than non-chim mate
    :param mate_tuple: output of return_mate_tuple
    :param chain: chain in ['+', '-']
    :return: distance, integer: >0 for "right" and < 0 for "left", None if chain not in ['+', '-'], =0 for overlap
    """
    chim_part1, chim_part2, mate2_dict, interval_intersection = mate_tuple
    mate1 = one_interval(dict_to_interval(chim_part1) | dict_to_interval(chim_part2))
    mate2 = one_interval(dict_to_interval(mate2_dict))
    right_dist = mate2[0].inf - mate1[0].sup
    left_dist = mate1[0].inf - mate2[0].sup
    sign = 0   # sign(configuration)
    if right_dist > 0:
        if chain == '+':
            sign = 1
        elif chain == '-':
            sign = -1
        else:
            return None
    elif left_dist > 0:
        if chain == '+':
            sign = -1
        elif chain == '-':
            sign = 1
        else:
            return None
    mate_dist = int(min(
        abs(right_dist),
        abs(left_dist),
    ))
    return mate_dist*sign


def reads_to_junctions(reads_df):
    """
    Takes dataframe with reads, groups by junctions
    :param reads_df: output of chim_input translated into dataframe
    :return: dataframe with junctions
    """
    index_list = ['chrom', 'chain', 'donor', 'acceptor']

    z = reads_df.groupby(index_list+['type']).size().reset_index(name='counts')
    zz = z.pivot_table(index=index_list,
                       columns=['type'], values=['counts'], fill_value=0)

    if not ('inside' in zz['counts'] or 'outside' in zz['counts']):
        if 'inside' in zz['counts']:
            zz['counts']['outside'] = 0
        if 'outside' in zz['counts']:
            zz['counts']['inside'] = 0
    zz = zz.sort_values([('counts', 'outside'), ('counts','inside')], ascending=False)

    xxd = reads_df.pivot_table(index=index_list,
                              values=['annot_donor'],
                              aggfunc='first')
    xxa = reads_df.pivot_table(index=index_list,
                              values=['annot_acceptor'],
                              aggfunc='first')
    chd = reads_df.pivot_table(index=index_list,
                              values=['chim_dist'],
                              aggfunc='first')
    yy = reads_df.pivot_table(index=index_list,
                              values=['id'],
                              aggfunc=lambda id: len(id.unique()))
    res = pd.concat([zz, xxd, xxa, chd, yy], axis=1)
    mi = res.columns
    ind = pd.Index([e[0] + '_'+ e[1] if e[0] == 'counts' else e for e in mi.tolist()])
    res.columns = ind

    return res.sort_values(['counts_outside','counts_inside'], ascending=False)

# Main


def main():
    # Arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="STAR output, Chimeric.out.junction \
                            OR list of such files")
    parser.add_argument("-s", "--sam", type=str,
                        help="Filtered STAR SAM output, with read_names same as in Chimeric.out.junction OR list")
    parser.add_argument("-o", "--output", type=str,
                        help="Output folder for results")
    parser.add_argument("-gz", "--gzip", type=str,
                        help="Option to create .json.gz")
    parser.add_argument("-l", "--list", type=str,
                        help="Enables list input mode. Options: input, sam, tag - MUST be lists")
    parser.add_argument("-gtf", "--gtf_annot", type=str,
                        default='/home/sunnymouse/Human_ref/hg19_exons.gtf',
                        help="Absolute path to annotation file")
    parser.add_argument("-t", "--tag", type=str,
                        default='ENCODE',
                        help="Tag name for grouping results, i.e. ENCODE id OR list of tags")
    args = parser.parse_args()

    # Exons GTF to junctions dict
    PTES_logger.info('Reading GTF...')
    gtf_exons_name = args.gtf_annot
    gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)
    PTES_logger.info('Reading GTF... done')

    # non-iterative
    make_dir(args.output)
    junc_dict = defaultdict(list)
    norm_junc_dict = defaultdict(lambda: defaultdict(list))

    if args.list:
        with open(args.input, 'r') as chim_names_file:
            chim_names_list = [x.strip('\n') for x in chim_names_file.readlines()]
        with open(args.sam, 'r') as sam_names_file:
            sam_names_list = [x.strip('\n') for x in sam_names_file.readlines()]
        with open(args.tag, 'r') as tag_names_file:
            tag_names_list = [x.strip('\n') for x in tag_names_file.readlines()]
        triads = zip(chim_names_list, sam_names_list, tag_names_list)
        PTES_logger.info('Enabled list mode')
        chim_reads_df_list = []
    else:
        triads = [(args.input, args.sam, args.tag)]

    for chim_name, sam_name, tag in triads:
        PTES_logger.info('Input file %s' % chim_name)

        # Reading filtered STAR non-chim output
        PTES_logger.info('Reading STAR non-chimeric output...')
        sam_dict = sam_input(sam_name=sam_name,
                            chim_name=chim_name,
                            norm_dict=norm_junc_dict)
        norm_read_names_list = []
        for k, v in norm_junc_dict.items():
            for read_name, read_dicts in v.items():
                norm_read_names_list.append({
                    'read_name' : read_name,
                    'chrom': k[0],
                    'chain': k[1],
                    'donor': k[2],
                    'acceptor': k[3],
                    'annot_donor': 1 if k[2] in gtf_donors[k[0]] else 0,
                    'annot_acceptor': 1 if k[3] in gtf_acceptors[k[0]] else 0,
                    'n_reads': len(v.items()),
                })
        try:
            norm_read_names = pd.DataFrame(norm_read_names_list)
        except KeyError:
            PTES_logger.warning('Creating norm split reads dataframe... empty dataframe')

        PTES_logger.info('Reading STAR non-chimeric output... done')

        # Reading filtered STAR output
        PTES_logger.info('Reading STAR chimeric output...')
        read_names_list = chim_input(chim_name=chim_name,
                                     gtf_donors=gtf_donors,
                                     gtf_acceptors=gtf_acceptors,
                                     sam_dict=sam_dict,
                                     junc_dict=junc_dict)
        PTES_logger.info('Reading STAR chimeric output... done')
        PTES_logger.info('Creating reads dataframes...')
        try:
            reads_df = pd.DataFrame(read_names_list)
            reads_df['id'] = tag
            if args.list:
                chim_reads_df_list.append(reads_df)
            else:
                all_reads_df = reads_df
        except KeyError:
            PTES_logger.warning('Creating reads dataframe... empty dataframe')

    if args.list:
        all_reads_df = pd.concat(chim_reads_df_list, sort=True).reset_index(drop=True)

    # Writing reads dataframe
    PTES_logger.info('Writing reads dataframes...')
    all_reads_df.to_csv(os.path.join(args.output, 'chim_reads.csv'), sep='\t')
    norm_read_names.to_csv(os.path.join(args.output, 'norm_split_reads.csv'), sep='\t')
    PTES_logger.info('Writing reads dataframes... done')

    # Writing junc_dict
    PTES_logger.info('Writing intervals to json files...')
    if args.gzip:
        PTES_logger.info('Output will be archived')
        with gzip.GzipFile(os.path.join(args.output, 'junc_dict.json.gz'), 'w') as junc_json, \
                gzip.GzipFile(os.path.join(args.output, 'norm_dict.json.gz'), 'w') as norm_json:
            junc_json.write(json.dumps({str(k): v for k, v in junc_dict.items()}).encode('utf-8'))
            norm_json.write(json.dumps({str(k1): v1 for k1, v1 in norm_junc_dict.items()}).encode('utf-8'))
    else:
        with open(os.path.join(args.output, 'junc_dict.json'), 'w') as junc_json, \
                open(os.path.join(args.output, 'norm_dict.json'), 'w') as norm_json:
            json.dump({str(k): v for k, v in junc_dict.items()}, junc_json, indent=2)
            json.dump({str(k1): v1 for k1, v1 in norm_junc_dict.items()}, norm_json, indent=2)

    PTES_logger.info('Writing intervals to json files... done')

    # Writing junctions dataframe
    PTES_logger.info('Creating junctions dataframe...')
    junctions_df = reads_to_junctions(reads_df=all_reads_df)
    junctions_df.to_csv(os.path.join(args.output, 'chim_junctions.csv'), sep='\t')
    PTES_logger.info('Creating junctions dataframe... done')


if __name__ == "__main__":
    main()
