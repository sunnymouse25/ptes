# Takes STAR Chimeric.out.junction output,
# converts circles into reads and junctions tables, writes intervals for junctions_to_bed into JSON file

# Imports

import argparse
from collections import defaultdict
import os
import json

import pandas as pd

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir, digit_code
from ptes.ptes import get_read_interval, star_line_dict, annot_junctions

# Functions


def reads_to_junctions(reads_df):
    index_list = ['chrom', 'chain', 'donor', 'acceptor']
    chd = reads_df.pivot_table(index=index_list,
                               values=['chim_dist', 'annot_donor', 'annot_acceptor', 'letters_ss'],
                               aggfunc='first')
    c = reads_df.pivot_table(index=index_list,
                             values=['read_name'],
                             aggfunc='count')
    c.columns=['counts']
    res = pd.concat([chd, c], axis=1)
    return res.sort_values(index_list)


def main():
    ### Arguments

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str,
                        help="STAR Chimeric.out.junction output")
    parser.add_argument("-o", "--output", type=str,
                        help="Path for subfolder with results")
    parser.add_argument("-f", "--filter", type=str,
                        help="Enable filtering for STAR output, creates .filtered file")
    parser.add_argument("-gtf", "--gtf_annot", type=str,
                        default='/home/sunnymouse/Human_ref/hg19_exons.gtf',
                        help="Absolute path to annotation file")
    parser.add_argument("-t", "--tag", type=str,
                        default='ENCODE',
                        help="Tag name for grouping results (prefix), i.e. ENCODE id")
    args = parser.parse_args()

    # Main
    make_dir(args.output)

    skipped = {'non-filtered': 0,    # different chromosomes and/or chains
               'chrM': 0,      # mapping to chrM
               'j_type-': 0,   # junction between the mates, -1 in STAR output
               'non-chim': 0}   # STAR counts very long (>1Mb) junctions as chimeric

    junc_dict = defaultdict(list)
    read_names_list = []

    if args.filter:
        PTES_logger.info('Filtering STAR output...')
        filtered_name = args.input.strip() + '.filtered'
        shell_call("cat %s | awk '$1 ==$4 && $3 ==$6' > %s" % (args.input, filtered_name))
        input_name = filtered_name
        PTES_logger.info('Filtering STAR output... done')
    else:
        input_name = args.input

    # Exons GTF to junctions dict
    PTES_logger.info('Reading GTF...')
    gtf_exons_name = args.gtf_annot
    gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=gtf_exons_name)
    PTES_logger.info('Reading GTF... done')

    annot_donors = 0
    annot_acceptors = 0

    PTES_logger.info('Input file: %s ' % input_name)

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
            junc_dict[(chrom, chain, line_dict['donor_ss'], line_dict['acceptor_ss'])]. append((chim_part1, chim_part2))

            annot_donor = 0
            annot_acceptor = 0
            if line_dict['donor_ss'] in gtf_donors[chrom]:
                annot_donor = 1
                annot_donors += 1
            if line_dict['acceptor_ss'] in gtf_acceptors[chrom]:
                annot_acceptor = 1
                annot_acceptors += 1

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
            }
            read_names_list.append(read_attrs)

    PTES_logger.info('Reading STAR output... done')
    PTES_logger.info('Processed: %i rows' % i)
    for key in skipped:
        PTES_logger.info('Skipped %s: %i rows' % (key, skipped[key]))
    PTES_logger.info('Converted successfully: %i rows' % len(read_names_list))
    PTES_logger.info('Annot donors: %i' % annot_donors)
    PTES_logger.info('Annot acceptors: %i' % annot_acceptors)
    PTES_logger.info('Creating reads dataframe...')
    try:
        reads_df = pd.DataFrame(read_names_list)
        reads_df['id'] = args.tag
        # Writing reads dataframe

        reads_df.to_csv(os.path.join(args.output, 'chim_reads.csv'), sep='\t')
        PTES_logger.info('Creating reads dataframe... done')

        # Writing junc_dict
        PTES_logger.info('Writing intervals to json...')
        with open(os.path.join(args.output, 'junc_dict.json'), 'w') as junc_json:
            json.dump({str(k): v for k, v in junc_dict.items()}, junc_json, indent=2)
        PTES_logger.info('Writing intervals to json... done')

        # Writing junctions dataframe
        PTES_logger.info('Creating junctions dataframe...')
        junctions_df = reads_to_junctions(reads_df=reads_df)
        junctions_df.to_csv(os.path.join(args.output, 'chim_junctions.csv'), sep='\t')
        PTES_logger.info('Creating junctions dataframe... done')
    except KeyError:
        PTES_logger.warning('Creating reads dataframe... empty dataframe')


if __name__ == "__main__":
    main()