# Takes .sam file
# For every split-read ('N' in CIGAR) makes JSON file with read_intervals for junctions_to_bed
# Input can be filtered by list of read_names in STAR Chimeric.out.junction output
# Also makes table with junctions, read_names and number of reads mapped to this junction

# Imports
from collections import defaultdict
import argparse
import os
import json
import gzip

import pandas as pd

from ptes.constants import PTES_logger
from ptes.lib.general import make_dir
from ptes.ptes import get_read_interval, parse_sam_row, annot_junctions


# Functions

def reads_to_junctions(reads_df, gtf_donors, gtf_acceptors):
    """
    Takes dataframe with reads, groups by junctions
    :param reads_df: reads dataframe
    :return: dataframe with junctions
    """
    index_list = ['chrom', 'chain', 'donor', 'acceptor']
    yy = reads_df.pivot_table(index=index_list,
                              values=['id'],
                              aggfunc=lambda id: len(id.unique()))
    yy['annot_donor'] = 0
    yy['annot_acceptor'] = 0
    for idx, row in yy.iterrows():
        row['annot_donor'] = 1 if idx[2] in gtf_donors[idx[0]] else 0
        row['annot_acceptor'] = 1 if idx[3] in gtf_acceptors[idx[0]] else 0

    yy.rename(index=str, columns={"id": "n_reads"})
    return yy.sort_values(index_list)


def main():
    # Arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="STAR output, Chimeric.out.junction.filtered \
                                OR list of such files")
    parser.add_argument("-s", "--sam", type=str,
                        help="Filtered STAR SAM output OR list")
    parser.add_argument("-o", "--output", type=str,
                        help="Output folder for results")
    parser.add_argument("-gz", "--gzip", type=str,
                        help="Option to create .json.gz")
    parser.add_argument("-l", "--list", type=str,
                        help="Enables list input mode. Options: sam, tag - MUST be lists")
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

    make_dir(args.output)
    norm_junc_dict = defaultdict(dict)
    norm_read_names_list = []

    if args.list:
        with open(args.sam, 'r') as sam_names_file:
            sam_names_list = [x.strip('\n') for x in sam_names_file.readlines()]
        if args.input:
            with open(args.input, 'r') as chim_names_file:
                chim_names_list = [x.strip('\n') for x in chim_names_file.readlines()]
        else:
            chim_names_list = [None]*len(sam_names_list)
        with open(args.tag, 'r') as tag_names_file:
            tag_names_list = [x.strip('\n') for x in tag_names_file.readlines()]
        triads = zip(chim_names_list, sam_names_list, tag_names_list)
        PTES_logger.info('Enabled list mode')
    else:
        triads = [(args.input, args.sam, args.tag)]

    for chim_name, sam_name, tag in triads:
        if chim_name:
            with open(chim_name, 'r') as chim_file:
                names_list = [x.strip('\n').split('\t')[9] for x in chim_file.readlines()]
                names_set = set(names_list)  # only names with chimeric output

        with open(sam_name, 'r') as sam_file:
            PTES_logger.info('Input file %s' % sam_name)
            for line in sam_file:
                if line.startswith('@'):
                    continue
                row = line.strip().split('\t')
                sam_attrs = None
                if len(row) > 1:
                    read_name = row[0]
                    if chim_name:
                        if read_name in names_set:
                            sam_attrs = parse_sam_row(row)
                    else:
                            sam_attrs = parse_sam_row(row)
                if sam_attrs:
                    if 'N' in sam_attrs['cigar']:  # read mapped with intron
                        read_dict = get_read_interval(cigar=sam_attrs['cigar'],
                                                      leftpos=sam_attrs['leftpos'],
                                                      output='dict')
                        if sam_attrs['chain'] == '+':
                            donor_ss = int(read_dict['N1'][0].inf - 1)  # counts first N as intron
                            acceptor_ss = int(read_dict['N1'][0].sup + 1)
                        elif sam_attrs['chain'] == '-':
                            donor_ss = int(read_dict['N1'][0].sup + 1)
                            acceptor_ss = int(read_dict['N1'][0].inf - 1)
                        norm_junc_dict[(sam_attrs['chrom'],
                                        sam_attrs['chain'],
                                        donor_ss,
                                        acceptor_ss)].update({read_name: tuple([read_dict])})
                        norm_read_names_list.append({
                            'read_name': read_name,
                            'chrom': sam_attrs['chrom'],
                            'chain': sam_attrs['chain'],
                            'donor': donor_ss,
                            'acceptor': acceptor_ss,
                            'id': tag})
    try:
        norm_read_df = pd.DataFrame(norm_read_names_list)
        norm_read_df = norm_read_df[
            ['read_name', 'chrom', 'chain',
             'donor', 'acceptor', 'id', ]
        ].sort_values(by=['chrom', 'chain', 'donor', 'acceptor']).reset_index(drop=True)
        PTES_logger.info('Writing reads dataframe...')
        norm_read_df.to_csv(os.path.join(args.output, 'norm_split_reads.csv'), sep='\t')
        PTES_logger.info('Writing reads dataframe... done')
    except KeyError:
        PTES_logger.warning('Creating norm split reads dataframe... empty dataframe')

    # Writing junc_dict
    PTES_logger.info('Writing intervals to json files...')
    if args.gzip:
        PTES_logger.info('Output will be archived')
        with gzip.GzipFile(os.path.join(args.output, 'norm_dict.json.gz'), 'w') as norm_json:
            norm_json.write(json.dumps({str(k1): v1 for k1, v1 in norm_junc_dict.items()}).encode('utf-8'))
    else:
        with open(os.path.join(args.output, 'norm_dict.json'), 'w') as norm_json:
            json.dump({str(k1): v1 for k1, v1 in norm_junc_dict.items()}, norm_json, indent=2)

    PTES_logger.info('Writing intervals to json files... done')

    # Writing junctions dataframe
    PTES_logger.info('Creating junctions dataframe...')
    junctions_df = reads_to_junctions(reads_df=norm_read_df, gtf_donors=gtf_donors, gtf_acceptors=gtf_acceptors)
    junctions_df.to_csv(os.path.join(args.output, 'norm_junctions.csv'), sep='\t')
    PTES_logger.info('Creating junctions dataframe... done')


if __name__ == "__main__":
    main()



