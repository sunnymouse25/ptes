# Takes segemehl output for simulated data and known junctions table
# Makes pivot table with true and false junctions
# Code from Jupyter notebook

### Arguments and imports
import argparse
import os
import errno
import multiprocessing as mp


import numpy as np
import pandas as pd


from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call, worker
from ptes.ptes import get_read_interval, one_interval, get_subseq, annot_junctions, \
    split_by_chimeric
from ptes.sim.sim import segemehl_to_intervals, intervals_to_junctions, get_detected, is_circ, map_dicts


parser = argparse.ArgumentParser()
parser.add_argument("-i","--input", type=str,
                    help="Segemehl output file; SAM file without header")
parser.add_argument("-o","--output", type=str,
                    help="Output folder for results")
parser.add_argument("-g","--genome", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/GRCh37.p13.genome.fa',
                    help="Absolute path to genome file")
parser.add_argument("-gtf","--gtf_annot", type=str,
                    default = '/uge_mnt/home/sunnymouse/Human_ref/hg19_exons.gtf',
                    help="Absolute path to genome file")
parser.add_argument("-t","--tag", type=str,
                    help="Tag name for grouping results, i.e. ENCODE id")
args = parser.parse_args()

# Functions


# Main
dirname = os.path.dirname(os.path.realpath(args.input))
if dirname == '':
    dirname = '.'
try:
    os.makedirs(args.output)
except OSError as exc:
    if exc.errno != errno.EEXIST:
        raise
    pass

path_to_file = args.output.rstrip('/')

col_names = ['read_name','strand','chrom','leftpos','cigar','NH','XI','XX','XY','XQ']
# col_nums = [0,1,2,3,5,14,19]
tag_list = ['NH','XI','XX','XY','XQ']
# Reading SAM input
PTES_logger.info('Reading SAM input...')
read_intervals, read_infos = segemehl_to_intervals(segemehl_outfile=args.input)

PTES_logger.info('Reading SAM input... done')

PTES_logger.info('Reading GTF...')

gtf_donors, gtf_acceptors = annot_junctions(gtf_exons_name=args.gtf_annot)
PTES_logger.info('Reading GTF... done')

PTES_logger.info('Creating junctions table...')
junc_list = intervals_to_junctions(read_intervals=read_intervals,
                                   read_infos=read_infos,
                                   gtf_donors=gtf_donors,
                                   gtf_acceptors=gtf_acceptors)   # from mapped read intervals to list of junctions
mapped_junc_df = pd.DataFrame(junc_list)
mapped_junc_df = mapped_junc_df[['read_name',
                                 'aln',
                                 'n_junctions',
                                 'chrom',
                                 'chain',
                                 'donor',
                                 'annot_donor',
                                 'acceptor',
                                 'annot_acceptor',
                                 'chimeric']].sort_values(by=['read_name','aln']).reset_index(drop=True)
gr = mapped_junc_df.groupby(['read_name','aln']).apply(lambda x: x.chimeric.any()).reset_index(name='chim_read')
del gr['aln']
mapped_junc_df = pd.merge(mapped_junc_df, gr, on='read_name').reset_index(drop=True)
mapped_junc_df['tag'] = args.tag   # ENCODE id for grouping results
mapped_junc_df.to_csv('%s/mapped_junc_df_segemehl.csv.gz' % path_to_file, sep = '\t', compression='gzip')

PTES_logger.info('Creating junctions table... done')

PTES_logger.info('Comparing mapped junctions with real junctions...')
params_list = []
junc_df_list = mp.Manager().list()
file_prefix = '/home/sunnymouse/projects/PTES/Single-read/read_simulation/remap_sjdb/'

for junc_df in pd.read_csv(file_prefix + 'junc_df_detected.csv', sep='\t', index_col=0, chunksize=1000):
    params_list.append({'function': get_detected,
                        'df': junc_df,
                        'mapped': mapped_junc_df,
                        'df_list': junc_df_list})
pool = mp.Pool(20)
pool.map(worker, params_list)
pool.close()
pool.join()

junc_df = pd.concat(list(junc_df_list), ignore_index=True)
junc_df['junction'] = junc_df.apply(is_circ, axis=1)
junc_df.to_csv('%s/new_junc_df_detected.csv' % path_to_file, sep = '\t')

true_map_dict, false_map_dict = map_dicts(real_junc_df=junc_df, mapped_junc_df=mapped_junc_df)
false_map_df = pd.DataFrame.from_dict(false_map_dict, orient='index')
false_map_df.columns = ['real', 'mapped']
false_map_df['correct'] = False
true_map_df = pd.DataFrame.from_dict(true_map_dict, orient='index')
true_map_df.columns = ['real', 'mapped']
true_map_df['correct'] = True
all_map_df = pd.concat([false_map_df, true_map_df])
all_map_df.to_csv('%s/all_map_df.csv' % path_to_file, sep = '\t')

info_df = pd.merge(junc_df, all_map_df, left_on='read_name',right_index=True).reset_index(drop=True)
del info_df['real']
x = info_df.groupby(['n_junctions','mapped','correct']).apply(lambda x: x.read_name.nunique()).reset_index(name='counts')
y = pd.pivot_table(x, index=['n_junctions','correct'], columns = ['mapped'],values=['counts'], fill_value=0, aggfunc=sum, margins=True)
y.to_csv('%s/pivot_table.csv' % path_to_file, sep='\t')
pivot_html_name = 'pivot_table.html'
init_file(pivot_html_name, folder=path_to_file)
writeln_to_file(y.to_html(), pivot_html_name, folder=path_to_file)

PTES_logger.info('Comparing junctions table with real junctions... done')