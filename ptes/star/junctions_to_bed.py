# Takes junctions table, junc_dict with read intervals as json file and conditions for creating BED files
# Conditions are used to send query to junctions table
# Creates 4 BED files with intervals with conditions specified
# For viewing in UCSC Genome Browser creates .coords.csv and .track files

# Imports
from collections import OrderedDict
import argparse
import os
import json
import gzip

import pandas as pd
from interval import interval

from ptes.constants import PTES_logger
from ptes.lib.general import shell_call, make_dir, digit_code
from ptes.ucsc.ucsc import get_track_list, to_bigbed, get_single_track


# Functions


def main():
    # Arguments
    '''
    args_s = ('-t ../tests/test_data/ptes/chim_reads_test.csv '
    '-j ../tests/test_data/ptes/junc_dict.json.gz '
    '-f ../tests/test_data/ptes/chim_junctions_test.csv '
    '-q letters_ss=="." '
    '-o ../tests/test_results/bed '
    '-p test '
    '-gz 1')
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--table", type=str,
                        help="DataFrame with data to create BED, required 4 columns are: chrom, strand, donor, acceptor")
    parser.add_argument("-sep", "--separator", type=str,
                        default='\t',
                        help="DataFrame separator, tab by default")
    parser.add_argument("-j", "--json", type=str,
                        help="JSON file with all read intervals as OrderedDicts")
    parser.add_argument("-gz", "--gzip", type=str,
                        help="Write anything to enable reading .json.gz")
    parser.add_argument("-q", "--query", type=str,
                        help="Conditions to filter junctions table, string as in pandas.DataFrame.query()")
    parser.add_argument("-f", "--filter", type=str,
                        help="DataFrame with (chrom, strand, donor, acceptor) to filter input table")
    parser.add_argument("-n", "--names", type=str,
                        nargs='+',
                        default=['chim1', 'chim2', 'mate2', ],
                        help="List of names for chim parts in BED name: [chim1, chim2, mate2]. \
                        Important: same order as parts in json values")
    parser.add_argument("-c", "--colors", type=str,
                        nargs='+',
                        default=['r', 'r', 'b', ],
                        help="List of colors for chim parts in BED name: [chim1, chim2, mate2]. \
                        Important: same order as parts in json values\
                        Colors: 'r', 'g', 'b' or in RGB code like '0,255,0'")
    parser.add_argument("-o", "--output", type=str,
                        default='bed',
                        help="Output folder for results, default is bed/")
    parser.add_argument("-p", "--prefix", type=str,
                        default='Output',
                        help="Prefix for all output files")
    parser.add_argument("-sort", "--sort", type=str,
                        help="Write anything to enable sorting BED files")
    parser.add_argument("-bb", "--bigbed", type=str,
                        help="Write anything to enable creating .bigBed files")
    args = parser.parse_args()
#    args = parser.parse_args(args_s.split(' '))

    PTES_logger.info('Reading input files...')
    make_dir(args.output)

    index_list = ['chrom', 'chain', 'donor', 'acceptor']
    input_df = pd.read_csv(args.table, sep=args.separator)

    for col in index_list:
        if col not in input_df.columns:
            PTES_logger.error('Input table does not contain required column %s ' % col)
            os._exit(1)

    if args.filter:   # filter by junctions
        filter_df = pd.read_csv(args.filter, sep=args.separator)
        for col in index_list:
            if col not in filter_df.columns:
                PTES_logger.error('Filter table does not contain required column %s ' % col)
                os._exit(1)
        cols_to_use = index_list + list(input_df.columns.difference(filter_df.columns))  # avoid repeating columns
        df_new = pd.merge(filter_df, input_df[cols_to_use], on=index_list, how='inner',)
    else:
        df_new = input_df

    if args.query:   # filter reads by conditions
        df_new = df_new.query(args.query)

    df_new.to_csv(os.path.join(args.output, 'df_filter.csv'), sep='\t')

    # Reading .json.gz file
    if args.gzip:
        with gzip.GzipFile(args.json, 'r') as fin:
            junc_dict = json.loads(fin.read().decode('utf-8'), object_pairs_hook=OrderedDict)
    else:
        junc_dict = json.load(open(args.json), object_pairs_hook=OrderedDict)

    len_read_dicts = len(junc_dict.values()[0].values())   # must be 3 for mate_inside/outside and 2 for circles
    if len(args.names) < len_read_dicts:
        PTES_logger.warning('List of names has less items than list of features in read_dicts!')
        part_names = [x[0]+str(x[1]) for x in list(zip(['part_']*len_read_dicts, range(1, len_read_dicts+1)))]
    else:
        part_names = args.names

    if len(args.colors) < len_read_dicts:
        PTES_logger.warning('List of colors has less items than list of features in read_dicts!')
        part_colors = ['r']*len_read_dicts
    else:
        part_colors = args.colors

    PTES_logger.info('Reading input files... done')
    PTES_logger.info('Creating BED files...')
    bed_name = '%s.bed' % args.prefix  # only track lines
    unique_bed_name = '%s.unique.bed' % args.prefix  # one representative read for unique junctions
    single_bed_name = '%s.single.bed' % args.prefix  # single line for one chimeric junction
    single_unique_bed_name = '%s.single.unique.bed' % args.prefix  # for unique junctions, single line for one junction
    code_name = '%s.codes.csv' % args.prefix   # table with codes for each read and descriptions
    coord_name = '%s.coords.csv' % args.prefix  # table with windows to paste into GB and descriptions
    info_name = '%s.track' % args.prefix  # file to submit to GB

    bed_list = []  # for outputting BED lines
    unique_dict = {}  # for outputting BED lines, unique chimeric junctions
    single_list = []  # for outputting BED lines, one row per one chimeric junction
    single_unique_list = []  # for outputting BED lines, one row per one unique chimeric junction
    coord_list = []  # for outputting coord lines
    code_list = []  # for outputting coord lines, one row per one read

    with open(os.path.join(args.output, info_name), 'w') as info_file:
        info_file.write('\n'.join(
        ['browser full knownGene ensGene cons100way wgEncodeRegMarkH3k27ac rmsk',
         'browser dense refSeqComposite pubs snp150Common wgEncodeRegDnaseClustered wgEncodeRegTfbsClusteredV3',
         'browser pack gtexGene',
         'track type=bigBed \
         name="%s" \
         description="bigBed" \
         visibility=2 \
         itemRgb="On" \
         bigDataUrl=https://github.com/sunnymouse25/ptes/blob/dev/research/bed/%s?raw=true' % (
             args.prefix,
             bed_name.replace('.bed', '.bb')
             )
            ]
           )
        )

    num = 0
    junctions = df_new.groupby(index_list)['read_name'].apply(list) # unique chimeric junctions
    for index, read_list in junctions.items():  # value is list of read_names
        chrom = index[0]  # index is (chrom, chain, donor_ss, acceptor_ss)
        chain = index[1]
        donor_ss = index[2]
        acceptor_ss = index[3]
        windows_min = []
        windows_max = []
        codes = []
        for read_name in read_list:  # for each read w. this junction
            num += 1
            code = digit_code(number=num)  # every unique number will be 6-digit
            codes.append(code)
            track_lists = []
            if not unique_dict.get(index, None):  # for each unique junction write the 1st read line(s)
                unique_dict[index] = []
                add_unique = True
            else:
                add_unique = False
            read_dict_list = junc_dict[str(index)][read_name]  # list of dicts: each dict is one track (i.e. chim_part)
            # Iterating over tracks
            for i, read_dict in enumerate(read_dict_list):
                for k, v in read_dict.items():
                    read_dict[k] = interval[v[0][0], v[0][1]]
                track_list = get_track_list(chrom=chrom,
                                            chain=chain,
                                            read_dict=read_dict,
                                            name='_'.join(map(str, [donor_ss, acceptor_ss, code, part_names[i]])),
                                            color=part_colors[i])
                track_lists.append(track_list)
            # Writing BED lines, collecting extremas for window size
            for track_list in track_lists:
                windows_min.append(int(track_list[1]))  # track_list[1] is chromStart, track_list[2] is chromEnd
                windows_max.append(int(track_list[2]))
                bed_line = '\t'.join(track_list)
                bed_list.append(bed_line)
                if add_unique:
                    unique_dict[index].append(bed_line)
            # Writing code line
            code_list.append({
                'chrom': chrom,
                'chain': chain,
                'donor': donor_ss,
                'acceptor': acceptor_ss,
                'read_name': read_name,
                'code': code
            })
            # Making BED file with one row for the pair of mates
            single_track = get_single_track(read_dict_list=read_dict_list,
                                            kwargs={'chrom': chrom,
                                                    'chain': chain,
                                                    'name': '_'.join(
                                                        map(str, [donor_ss, acceptor_ss, code])),
                                                    'color': '255,0,255'})  # for checking in GB that intervals are same
            single_list.append('\t'.join(single_track))
            if add_unique:
                single_unique_list.append('\t'.join(single_track))
        # Description for the junction into coords.csv
        window = (chrom,  # one window for junction
                  min(windows_min) - 200,
                  max(windows_max) + 200)
        coord_list.append({
                    'chrom': chrom,
                    'chain': chain,
                    'donor': donor_ss,
                    'acceptor': acceptor_ss,
                    'window': '%s:%i-%i' % window,
                    'codes': '-'.join(map(str,[codes[0], codes[-1]])),
        })

    PTES_logger.info('Creating BED files... done')

    PTES_logger.info('Writing BED files...')
    with open(os.path.join(args.output, bed_name), 'w') as bed_file, \
            open(os.path.join(args.output, unique_bed_name), 'w') as unique_bed_file, \
            open(os.path.join(args.output, single_bed_name), 'w') as single_bed_file, \
            open(os.path.join(args.output, single_unique_bed_name), 'w') as single_unique_bed_file, \
            open(os.path.join(args.output, coord_name), 'w') as coord_file, \
            open(os.path.join(args.output, code_name), 'w') as code_file:
        bed_file.write('\n'.join(bed_list))
        single_bed_file.write('\n'.join(single_list))
        single_unique_bed_file.write('\n'.join(single_unique_list))

        for unique_value in unique_dict.values():
            unique_bed_file.write('\n'.join(list(unique_value))+'\n')

    PTES_logger.info('Writing BED files... done')

    PTES_logger.info('Creating junctions dataframes...')
    coord_df = pd.DataFrame(coord_list)
    code_df = pd.DataFrame(code_list)
    coord_df.to_csv(os.path.join(args.output, coord_name), sep='\t')
    code_df.to_csv(os.path.join(args.output, code_name), sep='\t')

    PTES_logger.info('Creating junctions dataframes... done')

    if args.sort:
        PTES_logger.info('Sorting BED files...')
        for filename in [bed_name, unique_bed_name, single_bed_name, single_unique_bed_name]:
            shell_call('cat %s | sort -k1,1 -k2,2n  > %s.sorted' % (os.path.join(args.output, filename),
                                                                    os.path.join(args.output, filename),)                                                                                                                                                    )

        PTES_logger.info('Sorting BED files... done')

    if args.bigbed:   # will also sort files
        PTES_logger.info('Making bigBed...')
        for filename in [bed_name, unique_bed_name, single_bed_name, single_unique_bed_name]:
            to_bigbed(bed_name=filename, folder_name=args.output)
        PTES_logger.info('Making bigBed... done')


if __name__ == "__main__":
    main()