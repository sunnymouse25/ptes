from collections import defaultdict, OrderedDict
import os
import multiprocessing as mp

import numpy as np
import pandas as pd

from ptes.ptes import get_read_interval, get_interval_length, one_interval, sort_by_xq, get_junctions
from ptes.lib.general import init_file, writeln_to_file, shell_call, worker


def read_coord(intervals, infos, read_name, cigar, leftpos, NH, XI, XQ, chrom, flag):
    xq = int(XQ)
    if flag & 16 == 0:
        chain = '+'
    else:
        chain = '-'
    infos[read_name][XI].append(chrom)
    infos[read_name][XI].append(chain)
    read_interval = one_interval(
        get_read_interval(cigar, leftpos, output='interval'))  # local alignment is always bound by M
    intervals[read_name][XI].append((xq, read_interval))


def segemehl_to_intervals(segemehl_outfile):
    tag_list = ['NH', 'XI', 'XX', 'XY', 'XQ']
    read_intervals = defaultdict(lambda: defaultdict(list))  # mapped intervals
    read_infos = defaultdict(lambda: defaultdict(list))  # mapped chrom(s) and chain(s)
    with open(segemehl_outfile, 'r') as df_segemehl:
        for line in df_segemehl:
            row = line.strip().split('\t')
            read_name = row[0]
            flag = int(row[1])
            chrom = row[2]
            leftpos = row[3]
            cigar = row[5]
            sam_attrs = {'read_name': read_name,
                         'flag': flag,
                         'chrom': chrom,
                         'leftpos': leftpos,
                         'cigar': cigar}
            tags = dict.fromkeys(tag_list, None)
            for elm in row[13:]:
                for tag in tag_list:
                    if tag in elm:
                        tags[tag] = elm
            sam_attrs.update(tags)
            sam_attrs['NH'] = sam_attrs['NH'].strip('NH:i:')
            if sam_attrs['NH'] != '1':
                continue   # unique mappings only
            sam_attrs['XI'] = sam_attrs['XI'].strip('XI:i:')
            if sam_attrs['XQ']:
                sam_attrs['XQ'] = sam_attrs['XQ'].strip('XQ:i:')
                read_coord(read_intervals, read_infos, **sam_attrs)
    return read_intervals, read_infos


def intervals_to_junctions(read_intervals, read_infos, gtf_donors, gtf_acceptors):
    junc_list = []  # from mapped read intervals to the list of junctions
    for key in read_intervals:  # key is read_name
        xi = 0
        while read_intervals[key].get(str(xi), None):  # xi - number of current read alignment
            tuples = read_intervals[key][str(xi)]  # list of intervals of current read alignment
            infos = read_infos[key][str(xi)]
            chroms = set(infos[::2])
            chains = set(infos[1::2])
            if len(chroms) == 1 and len(chains) == 1:  # read must be mapped to the same chrom and chain
                chrom = chroms.pop()
                if chrom.startswith('chr'):  # skip contigs
                    chain = chains.pop()
                    values = sort_by_xq(tuples=tuples, chain=chain)
                    xi += 1
                    if len(values) > 0:  # should be every read, just to be sure
                        attrs = {'read_name': key,
                                'aln': xi-1}
                        for junction in get_junctions(chrom=chrom,
                                                      chain=chain,
                                                      values=values,
                                                      gtf_donors=gtf_donors,
                                                      gtf_acceptors=gtf_acceptors):
                            attrs.update(junction)
                            junc_list.append(attrs)
    return junc_list


def detect_ss(row, mapped):
    if row['n_junctions'] > 0:
        if mapped.query('read_name == "%s" & donor == "%s" & acceptor == "%s" '
                        % (row['read_name'], row['donor'], row['acceptor'])).shape[0] > 0:
            return True
        elif mapped.query('read_name == "%s" & donor == "%s" & acceptor == "%s" '
                          % (row['read_name'], row['acceptor'], row['donor'])).shape[0] > 0:
            return True
        else:
            return False
    else:
        temp = mapped.query('read_name == "%s" & chrom == "%s"' % (row['read_name'], row['chrom']))
        if temp.shape[0] > 0:
            if temp.isnull().values.any():
                return True
            else:
                return False
        else:
            return False


def get_detected(df, mapped, df_list, **kwargs):
    df.donor = df.donor.dropna().apply(lambda x: str(int(x)))
    df.acceptor = df.acceptor.dropna().apply(lambda x: str(int(x)))
    del df['detected']
    df['detected'] = df.apply(detect_ss, args=(mapped,),axis=1)
    junc_df_list = df_list
    junc_df_list.append(df)

def is_circ(row):
    """
    for table with real junction written in 'exons' column
    :param row: should have 'exons' column separated by ','
    :return: type of junction: ordinary or circular
    """
    if isinstance(row['exons'], str):
        exons = row['exons'].partition(',')
        ex1 = int(exons[0])
        ex2 = int(exons[2])
        if ex1<ex2:
            return 'ordinary'
        elif ex1>ex2:
            return 'circular'
    else:
        return np.nan


def get_correct_reads(row, mapped, real_junc_df, true_map_dict, false_map_dict):
    """
    Makes reads with all their junctions mapped correctly (true_map_dict)
    and all other reads (false_map_dict)
    No paralleling
    :param row: row of df with real junctions and 'detected' column
    :param mapped: df with mapped junctions
    :param real_junc_df: df with real junctions
    :param true_map_dict: reads with all their junctions mapped correctly
    :param false_map_dict: all other reads
    :return: nothing, just updates 2 dictionaries
    """
    read_name = row['read_name']
    real_junc = row['n_junctions']
    mapped_junc_df = mapped   # global name of the df with mapped junctions
    if not row.detected:   # if any False than read goes to false_map_dict
        if not false_map_dict.get(read_name, None):              # if not there already is this read_name
            query = mapped_junc_df.query('read_name == "%s"' % read_name)
            if query.shape[0] > 0:
                mapped_junc = query.n_junctions.iloc[0]
            else:
                mapped_junc = 0
            false_map_dict[read_name] = [real_junc, mapped_junc]
    else:   # we should check if all other junctions are true
        if real_junc_df.query('read_name == "%s" & n_junctions == "%s"' % (read_name, str(real_junc))).detected.all():
            if not true_map_dict.get(read_name, None):
                query = mapped_junc_df.query('read_name == "%s"' % read_name)
                if query.shape[0] > 0:
                    mapped_junc = query.n_junctions.iloc[0]
                else:
                    mapped_junc = 0
                true_map_dict[read_name] = [real_junc, mapped_junc]   # just to be sure


def map_dicts(real_junc_df, mapped_junc_df):
    """
    Runs get_correct_reads function
    :param real_junc_df:
    :param mapped_junc_df:
    :return: reads with all their junctions mapped correctly (true_map_dict)
    and all other reads (false_map_dict)
    """
    true_map_dict = mp.Manager().dict()  # for reads with all their junctions mapped correctly
    false_map_dict = mp.Manager().dict()  # for all other reads
    real_junc_df.apply(get_correct_reads, args=(mapped_junc_df,real_junc_df, true_map_dict, false_map_dict), axis=1)
    return true_map_dict, false_map_dict

