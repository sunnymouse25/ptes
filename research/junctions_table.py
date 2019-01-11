# Makes actual junctions table using .aln and annotation for simulated reads

### Arguments and imports
import pandas as pd
import numpy as np
import subprocess, sys
import copy
from collections import defaultdict, OrderedDict

import numpy as np
import pandas as pd

from ptes.constants import PTES_logger
from ptes.lib.general import init_file, writeln_to_file, shell_call, worker
from ptes.ptes import get_read_interval, one_interval, get_subseq, annot_junctions, \
    split_by_chimeric
from ptes.sim.sim import segemehl_to_intervals, intervals_to_junctions, get_detected, is_circ, map_dicts


### FUNCTIONS


### MAIN

# Reading GTF
chrom_dict = {}
chain_dict = {}
start_dict = defaultdict(dict)
end_dict = defaultdict(dict)
with open('tr_names.gtf', 'r') as gtf_file:
    for line in gtf_file:
        line_list = line.strip().split('\t')
        chrom = line_list[0]
        strt = int(line_list[3])
        end = int(line_list[4])
        chain = line_list[6]
        attr = line_list[8]
        attr_list = attr.split(';')
        for attr in attr_list:
            if 'transcript_id' in attr:
                tr_name = attr.lstrip('transcript_id "').rstrip('"')
            elif 'exon_number' in attr:
                exon_number = attr.lstrip('exon_number ').strip()
        chain_dict[tr_name] = chain
        chrom_dict[tr_name] = chrom
        start_dict[tr_name][exon_number] = strt
        end_dict[tr_name][exon_number] = end

# Reading read table

sim_reads = pd.read_csv('sim_reads_results.csv', index_col=False, sep='\t')
sim_reads = sim_reads.set_index(["0"])

# Making the actual table - from annotation

x = 0  # will be index for this dataframe
for index, row in sim_reads.iterrows():
    exon_list = row['3'].split(',')
    n_exons = len(exon_list)
    n_junc = n_exons - 1
    tr_name = index.partition('_')[0]
    chrom = chrom_dict[tr_name]
    chain = chain_dict[tr_name]
    if chain == '+':
        for i, exon in enumerate(exon_list):
            if i == 0:
                prev_exon = exon
                try:
                    prev_donor = end_dict[tr_name][exon] + 1
                except KeyError:
                    print tr_name, exon
                    print type(exon)
            else:
                acc = start_dict[tr_name][exon] - 1
                junc_df.loc[x] = [index, n_junc, ','.join([prev_exon, exon]), chrom, chain, prev_donor, acc]
                prev_exon = exon
                prev_donor = end_dict[tr_name][exon] + 1
    if chain == '-':
        for i, exon in enumerate(exon_list):
            if i == 0:
                prev_exon = exon
                try:
                    prev_donor = start_dict[tr_name][exon] - 1
                except KeyError:
                    print tr_name, exon
                    print type(exon)
            else:
                acc = end_dict[tr_name][exon] + 1
                junc_df.loc[x] = [index, n_junc, ','.join([prev_exon, exon]), chrom, chain, prev_donor, acc]
                prev_exon = exon
                prev_donor = start_dict[tr_name][exon] - 1
    x += 1
